function [is_blackout,relay_outages,MW_lost,p_out,busessep,movie_data,n_msg] = dcsimsep(ps,br_outages,bus_outages,opt)  
% usage: [is_blackout,relay_outages,MW_lost,p_out,busessep,movie_data,n_msg] = dcsimsep(ps,br_outages,bus_outages,opt)  
% is_blackout indicates whether a large separation occurs
% branches_lost gives the set of dependant outages that occur due to relay actions
%  this is a ? by 2 matrix, with t in the first column, br no. in the second
% bus_outages gives the bus indices associated with bus failures
% MW_lost indicates how much load was lost as a result of small separations
% p_out is proportion of buses separated  
% busessep is a list of the buses that separated  
% n_msg is an nx1 vector of the maximum number of messages per iteration
% for each agent during the simulation (in case of distributed control)

% check the inputs
if nargin<3
    bus_outages = [];
end
if nargin<4
    opt = psoptions;
end

% init the outputs
is_blackout = 0;
MW_lost = struct('rebalance',0,'control',0);
imbalance = 0;
% rebalance: MW_lost after islanding in rebalance.m
% control: MW_lost due to load shedding in emergency control (central and
% distributed)
verbose = opt.verbose;

p_out=0;  
busessep=[];  

% Initialize
t = 0;
C = psconstants;
Pd0 = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
%Pd = Pd0;
Pd0_sum = sum(Pd0);
relay_outages = zeros(0,2);
ps.relay = relay_settings(ps,false,true);
% some constants
t_max = opt.sim.t_max; % time limit for the simulation
EPS = 1e-4;
dt_max = opt.sim.dt; % maximum simulation time

% set up agents if this is distributed
if strcmp(opt.sim.control_method,'distributed_control')
    ps_agents = set_up_agents(ps,opt);
    if verbose, fprintf('setting up agents...\n'); end
end

% Grab some useful data
nbus = size(ps.bus,1);
%m = size(ps.branch,1);
F = ps.bus_i(ps.branch(:,1));
T = ps.bus_i(ps.branch(:,2));
G = ps.bus_i(ps.gen(:,1));   % gen index
ge_status = ps.gen(:,C.ge.status);
Pg_max = ps.gen(:,C.ge.Pmax).*ge_status + EPS;
Pg_min = ps.gen(:,C.ge.Pmin).*ge_status - EPS;
Pg   = ps.gen(:,C.ge.Pg).*ge_status;
Pg0_sum = sum(Pg);
D = ps.bus_i(ps.shunt(:,1)); % load index
%NO_SEP = 0;
BIG_SEP = 2;
%SMALL_SEP = 1;
% set the power plant ramp rates
ramp_rate = ps.gen(:,C.ge.ramp_rate_up)/60; % ramp rate in MW/second
if all(ramp_rate==0)
    ramp_rate_MW_per_min = max(1,Pg_max*.05); % assume that all plants can ramp at 5% per minute. 
                                            % for a 100 MW plant, this
                                            % would be 5 MW/min. Reasonable
    ramp_rate = ramp_rate_MW_per_min/60;
end
ramp_dt = opt.sim.fast_ramp_mins*60;
ramp_limits = ramp_rate * ramp_dt;
ps.gen(:,C.ge.ramp_rate_down) = ramp_limits;

% Print the time
if verbose
    fprintf('------- t = 0.00 ----------\n');
end
% Step 1. rebalance and run the DCPF
if ~isfield(ps,'bus_i')
    ps = updateps(ps);
end
br_st = ps.branch(:,C.br.status)~=0;
% check to make sure that the base case is load balanced
if opt.debug && abs(Pd0_sum - Pg0_sum)>EPS
    error('The base case power system is not load balanced');
end
[sub_grids,n_sub_old] = find_subgraphs(ps.bus(:,1),ps.branch(br_st,1:2));
if n_sub_old>1
    error('The base case has more than one island');
end
% Find the ramp rate
ramp_rate( ~ge_status ) = 0; % plants that are shut down cannot ramp
% Check the mismatch
mis = total_P_mismatch(ps);
if abs(mis)>EPS, error('Base case has mismatch'); end
% Calculate the power flow
ps = dcpf(ps,[]); % this one should not need to do any rebalance, just line flow calcs
% Get the power flow
% flow = ps.branch(:,C.br.Pf);
% Record the movie data if requested

if nargout>5
    movie_data = record_movie_data([],t,ps,sub_grids);
end
% Error check
Pg = ps.gen(:,C.ge.Pg);
if opt.debug && any( Pg<Pg_min | Pg>Pg_max )
    error('Pg is out of bounds');
end

% Step 2. Apply exogenous outages
t = 1;
if verbose
    fprintf('------- t = %.3f ----------\n',t);
    fprintf('Exogenous events:\n');
end
% Apply the branch outages
if ~isempty(br_outages)
    ps.branch(br_outages,C.br.status) = 0;
    if verbose
        fprintf(' Removed branch %d\n',br_outages);
    end
end
% Apply the bus outages
if ~isempty(bus_outages)
    for i=1:length(bus_outages)
        bus_no = bus_outages(i);
        bus_ix = ps.bus_i(bus_no);
        if opt.debug && isempty(bus_ix) || bus_ix<=0 || bus_ix>nbus
            error('%d is not a valid bus number',bus_no);
        end
        br_set = (F==bus_ix) | (T==bus_ix);
        ps.branch(br_set,C.br.status) = 0;
        % trip gens and shunts at this bus
        ps.gen  (G==bus_ix,C.ge.status) = 0;
        ps.shunt(D==bus_ix,C.sh.status) = 0;
        if verbose
            fprintf(' Removed bus %d\n',bus_no);
        end
    end
end

% Record the movie data if requested
if nargout>5
    movie_data = record_movie_data(movie_data,t,ps,sub_grids);
end

% keep track of branches that have tripped
tripped_branches = [];

% Begin the main while loop for DCSIMSEP
it_no = 1;
t_prev_control = -dt_max; % the time that a previous control action was done.
                          % Set negative so that the first control action will
                          % be done no matter its time.
while true
    % Step 3. Find sub-grids in the network and check for major separation
    [sep,sub_grids,n_sub,p_out,busessep] = check_separation(ps,opt.sim.stop_threshold,verbose);
        
    % Step 4. rebalance & run the power flow
    %  if there are new islands, rebalance the generators
    if n_sub>n_sub_old
        ramp_dt = opt.sim.fast_ramp_mins*60; % do a minute of ramping before gen tripping/load shedding 
        max_ramp = ramp_rate*ramp_dt;
        Pd_old = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
        [Pg,ge_status,d_factor] = rebalance(ps,sub_grids,max_ramp,opt);
        % Error check:
        Pg_max = ps.gen(:,C.ge.Pmax).*ge_status + EPS;
        Pg_min = ps.gen(:,C.ge.Pmin).*ge_status - EPS;
%         if opt.debug && any( Pg<Pg_min | Pg>Pg_max ), error('Pg is out of bounds'); end
        if any( round(Pg)<round(Pg_min-EPS) | round(Pg)>round(Pg_max+EPS) ), error('Pg is out of bounds'); end
        % Implement the changes to load and generation
        ps.shunt(:,C.sh.factor) = d_factor;
        ps.gen(:,C.ge.status) = ge_status;
        ps.gen(:,C.ge.P) = Pg;
        ramp_rate(~ge_status) = 0; % make sure that failed generators don't ramp
        Pd_new = sum(ps.shunt(:,C.sh.P).*d_factor);
        MW_lost.rebalance = MW_lost.rebalance + Pd_old - Pd_new;
    end
    n_sub_old = n_sub;
    % run the power flow and record the flow
    ps = dcpf(ps,sub_grids);
    if opt.debug
        % Check that
        ge_status = ps.gen(:,C.ge.status);
        Pg_max = ps.gen(:,C.ge.Pmax).*ge_status + EPS;
        Pg_min = ps.gen(:,C.ge.Pmin).*ge_status - EPS;
        Pg = ps.gen(:,C.ge.Pg);
        if any( round(Pg)<round(Pg_min-EPS) | round(Pg)>round(Pg_max+EPS) ), 
            keyboard
            error('Pg is out of bounds');
        end
    end
    % Extract and record the flows
    %flow  = ps.branch(:,C.br.Pf);
    %br_st = ps.branch(:,C.br.status);
    if nargout>5
        movie_data = record_movie_data(movie_data,t,ps,sub_grids);
    end

    % Step 4a. Take control actions if needed.
    % Make sure the previous control action was done at least dt_max seconds ago
    if ~strcmp(opt.sim.control_method, 'none') && ...
            t - t_prev_control >= dt_max
        switch opt.sim.control_method
            % Compute and implement emergency control actions
            %  Note that this code also interfaces with a python comm. model,
            %  if requested in the options structure
            case 'emergency_control'
                [ps, this_MW_lost, mis] = central_control(ps,sub_grids,ramp_rate,opt,'dc');
            case 'distributed_control' % distributed emergency control 
                [ps_agents, ps, this_MW_lost, mis] = distributed_control(ps_agents,ps,sub_grids,ramp_rate,it_no,opt,'dc');
            otherwise
                error('Undefined control method.')
        end
        MW_lost.control = MW_lost.control + this_MW_lost.control;
        imbalance = imbalance + mis;
        % update previous control time
        t_prev_control = t;
    end
    
    % Step 5. Update relays
    [ps.relay,br_out_new,dt,n_over] = update_relays(ps,verbose,dt_max);
    multi_trip = ismember(br_out_new,tripped_branches);
    if any(multi_trip)
        error('Lines were tripped multiple times');
    end
    tripped_branches = union(tripped_branches,br_out_new);
    
    % Step 6. Check for any remaining overload potential, decide if we
    % should stop the simulation
    if n_over == 0 % previously: dt==Inf (changed because dt is now reduced
                   % to make sure control actions are implemented every minute)
        if verbose
            fprintf(' There are no overloads in the system. Quitting...\n');
        end
        break
    elseif n_over > 0 && verbose
        fprintf(' There are %d overloads in the system.\n',n_over);
    end
    
    % If we want to stop when the network is divided into subnetworks, do
    % this:
    if opt.sim.stop_on_sep
        if sep==BIG_SEP
            is_blackout = 1;
            if verbose
                fprintf('-------------- t = %.3f ----------------\n',t);
                fprintf('----------- Major separation -----------\n');
            end
            break
        end
    else % If we wanted to stop after a certain amount of load shedding, do this:
        Pd_sum = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
        load_remaining_fraction = Pd_sum/Pd0_sum;
        if verbose
            %fprintf('------------- t = %.3f ---------------\n',t);
            fprintf('-------- %.1f%% of load remains --------\n',load_remaining_fraction*100);
        end
        if load_remaining_fraction < opt.sim.stop_threshold && ~is_blackout
            is_blackout = 1;
            if verbose
                fprintf('----------- Blackout occurred ----------\n');
            end
        end
        if load_remaining_fraction<opt.sim.stop_threshold
            break
        end
    end
    
    % advance/print the time
    t = t + dt;
    if t > t_max % stop if the next outage happens after t_max
        t = t_max;
        break
    end
    if verbose
        fprintf('------- t = %.3f ----------\n',t);
    end
    
    % Step 7. Trip overloaded branches
    ps.branch(br_out_new,C.br.status) = 0;
    ps.branch(br_out_new,C.br.Imag) = 0;
    ps.branch(br_out_new,C.br.Pf) = 0;
    ps.branch(br_out_new,C.br.Pt) = 0;
    % record which branches were lost
    for i = 1:length(br_out_new)
        br = br_out_new(i);
        relay_outages = cat(1,relay_outages,[t br]);
    end
    
    % Record movie data
    if nargout>5
        movie_data = record_movie_data(movie_data,t,ps,sub_grids);
    end
    
    % print something
    if verbose && ~isempty(br_out_new)
        fprintf(' Branch %d tripped on overcurrent\n',br_out_new);
    end
    
    % Increment the counter and return to step 3.
    it_no = it_no + 1;
end
max_ramp = ramp_rate * opt.sim.fast_ramp_mins*60;
% do a final rebalance just to make sure
[Pg,ge_status,d_factor] = rebalance(ps,sub_grids,max_ramp,opt);
% Error check
% if opt.debug
    Pg_max = ps.gen(:,C.ge.Pmax).*ge_status + EPS;
    Pg_min = ps.gen(:,C.ge.Pmin).*ge_status - EPS;
%     if any( Pg<Pg_min | Pg>Pg_max ), error('Pg is out of bounds'); end
    if any( round(Pg)<round(Pg_min-EPS) | round(Pg)>round(Pg_max+EPS) ), error('Pg is out of bounds'); end
% end
% Implement
ps.shunt(:,C.sh.factor) = d_factor;
ps.gen(:,C.ge.status) = ge_status;
ps.gen(:,C.ge.P) = Pg;
% Compute the amount of load lost
Pd = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
total_MW_lost = Pd0_sum - sum(Pd);
if abs(total_MW_lost - (MW_lost.control + MW_lost.rebalance)) > EPS
    error('Something is wrong.')
end

if strcmp(opt.sim.control_method,'distributed_control')
    n_msg = zeros(nbus,1);
    % Compute the maximum number of messages per all iterations for each agent
    for i = 1:nbus
        n_msg(i) = max(ps_agents(i).messages);
    end
else
    n_msg = 0;
end

% Print something
if verbose
    n_overloads = sum(abs(ps.branch(:,C.br.Pf))>ps.branch(:,C.br.rateB)+EPS);
    fprintf('-------------- t = %7.3f -----------------\n',t);
    fprintf(' Simulation complete\n');
    fprintf('  %d emergency (rateB) overloads remain\n',n_overloads);
    fprintf('  %d endogenous relay outages\n',size(relay_outages,1));
    fprintf('  %g MW load lost (%.1f%%)\n',total_MW_lost,total_MW_lost/Pd0_sum*100);
    fprintf('  %g MW (%.1f%%) in rebalance,  %g MW (%.1f%%) in control\n', ...
        MW_lost.rebalance, MW_lost.rebalance/Pd0_sum*100, ...
        MW_lost.control, MW_lost.control/Pd0_sum*100);
    if imbalance > EPS
        fprintf('  %g MW (%.1f%%) imbalance occured due to control.',...
            imbalance, imbalance/Pd0_sum*100);
    end
    fprintf('--------------------------------------------\n');
end

% Record final movie data
if nargout>5
    movie_data = record_movie_data(movie_data,t,ps,sub_grids);
end

% Extra functions
function movie_data = record_movie_data(movie_data,t,ps,sub_grids)
% used to record data needed for movie making
C = psconstants;

% get data
Pg = ps.gen(:,C.ge.Pg).*ps.gen(:,C.ge.status);
Pd = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
br_st = ps.branch(:,C.br.status)~=0;
flow = ps.branch(:,C.br.Pf);

% Record data
if isempty(movie_data)
    movie_data.t = t;
    movie_data.flows = flow;
    movie_data.br_status = br_st;
    movie_data.sub_grids = sub_grids;
    movie_data.Pg = Pg;
    movie_data.Pd = Pd;
    movie_data.Pd0 = Pd;
    movie_data.Pg0 = Pg;
    movie_data.MW_lost = 0;
else
    movie_data.t = [movie_data.t t];
    movie_data.flows = [movie_data.flows flow];
    movie_data.br_status = [movie_data.br_status br_st];
    movie_data.sub_grids = [movie_data.sub_grids sub_grids];
    movie_data.Pg = [movie_data.Pg Pg];
    movie_data.Pd = [movie_data.Pd Pd];
    MW_lost = sum(movie_data.Pd0) - sum(Pd);
    movie_data.MW_lost = [movie_data.MW_lost MW_lost];
end
