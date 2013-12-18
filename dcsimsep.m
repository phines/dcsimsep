function [is_blackout,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,bus_outages,opt)  
% usage: [is_blackout,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,bus_outages,opt)  
% is_blackout indicates whether a large separation occurs
% branches_lost gives the set of dependant outages that occur due to relay actions
%  this is a ? by 2 matrix, with t in the first column, br no. in the second
% bus_outages gives the bus indices associated with bus failures
% MW_lost indicates how much load was lost as a result of small separations
% p_out is proportion of buses separated  
% busessep is a list of the buses that separated  

% check the inputs
if nargin<3
    bus_outages = [];
end
if nargin<4
    opt = psoptions;
end

% init the outputs
C = psconstants;
is_blackout = 0;

p_out=0;  
busessep=[];  

Pd0 = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
Pd0_sum = sum(Pd0);
relay_outages = zeros(0,2);
ps.relay = relay_settings(ps,false,true);
% some constants
dt_max = opt.sim.dt_max_default;
t_max = 60*30;
EPS = 1e-4;
% 
pid = feature('getpid');

% comm status file name for emergency control
comm_status_file = sprintf('/tmp/comm_status_%d.csv',pid);
grid_status_file = sprintf('/tmp/grid_status_%d.csv',pid);

% grab useful data
C = psconstants;
n = size(ps.bus,1);
m = size(ps.branch,1);
F = ps.bus_i(ps.branch(:,1));
T = ps.bus_i(ps.branch(:,2));
flow_max = ps.branch(:,C.br.rateB);
G = ps.bus_i(ps.gen(:,1));   % gen index
Pmax = ps.gen(:,C.ge.Pmax);
ge_status = ps.gen(:,C.ge.status);
Pg   = ps.gen(:,C.ge.Pg).*ge_status;
Pg0_sum = sum(Pg);
Pg_max = ps.gen(:,C.ge.Pmax).*ge_status;
Pg_min = ps.gen(:,C.ge.Pmin).*ge_status;
D = ps.bus_i(ps.shunt(:,1)); % load index
%NO_SEP = 0;
BIG_SEP = 2;
%SMALL_SEP = 1;
% set the power plant ramp rates
ramp_rate = ps.gen(:,C.ge.ramp_rate_up)/60; % ramp rate in MW/second
if all(ramp_rate==0)
    ramp_rate_MW_per_min = max(1,Pmax*.05); % assume that all plants can ramp at 5% per minute. 
                                            % for a 100 MW plant, this
                                            % would be 5 MW/min. Reasonable
    ramp_rate = ramp_rate_MW_per_min/60;
end
% emergency ramp rate is faster
ramp_rate_emergency = ramp_rate;
% print the time:
if opt.verbose
    fprintf('------- t = 0.00 ----------\n');
end
% Step 1. redispatch and run the DCPF
ps = updateps(ps);
br_st = ps.branch(:,C.br.status)~=0;
% check to make sure that the base case is load balanced
if abs(Pd0_sum - Pg0_sum)>EPS
    error('The base case power system is not load balanced');
end
[sub_grids,n_sub_old] = findSubGraphs(ps.bus(:,1),ps.branch(br_st,1:2));
if n_sub_old>1
    error('The base case has more than one island');
end
%ramp_rate( (Pg<EPS) | (~ge_status) ) = 0; % plants that are shut down cannot ramp
ramp_rate( ~ge_status ) = 0;
% calculate the power flow
ps = dcpf(ps,[],false); % this one should not need to do any redispatch, just line flow calcs
flow = ps.branch(:,C.br.Pf);
if nargout>5
    flows = flow;
end

% error check
if any(ps.gen(:,C.ge.Pg)>(ps.gen(:,C.ge.Pmax)+EPS)), disp(br_outages); error('Pg is out of bounds'); end

% Step 2. Apply exogenous outages
t = 1;
if opt.verbose
    fprintf('------- t = %.3f ----------\n',t);
    fprintf('Exogenous events:\n');
end
% branch outages
if ~isempty(br_outages)
    ps.branch(br_outages,C.br.status) = 0;
    if opt.verbose
        fprintf(' Removed branch %d\n',br_outages);
    end
end
% bus outages
if ~isempty(bus_outages)
    for i=1:length(bus_outages)
        bus_no = bus_outages(i);
        bus_ix = ps.bus_i(bus_no);
        if isempty(bus_ix) || bus_ix==0
            error('%d is not a valid bus number',bus_no);
        end
        br_set = (F==bus_ix) | (T==bus_ix);
        ps.branch(br_set,C.br.status) = 0;
        % trip gens and shunts at this bus
        ps.gen  (G==bus_ix,C.ge.status) = 0;
        ps.shunt(D==bus_ix,C.sh.status) = 0;
        if opt.verbose
            fprintf(' Removed bus %d\n',bus_no);
        end
    end
end

dt = 10.00; % This initial time step sets the quantity of initial gen ramping.
while t < t_max
    % Step 3. Find sub-grids in the network and check for major separation
    [sep,sub_grids,n_sub,p_out,busessep] = check_separation(ps,opt.sim.stop_threshold,opt.verbose);
    
    % Step 4. redispatch & run the power flow
    % if there are new islands, redispatch the generators
    if n_sub>n_sub_old
        ramp_dt = max(dt,opt.sim.fast_ramp_mins*60); % the amount of generator ramping time to allow
        max_ramp = ramp_rate*ramp_dt;
        [Pg,ge_status,d_factor] = redispatch(ps,sub_grids,max_ramp,opt.verbose);
        % error check:
        if any(Pg>ps.gen(:,C.ge.Pmax)+EPS), error('Pg is out of bounds'); end
        % implement the changes to load and generation
        ps.shunt(:,C.sh.factor) = d_factor;
        ps.gen(:,C.ge.status) = ge_status;
        ps.gen(:,C.ge.P) = Pg;
        ramp_rate(~ge_status)=0; % make sure that failed generators don't ramp
    end
    n_sub_old = n_sub;
    % run the power flow and record the flow
    ps = dcpf(ps,sub_grids,true);
    if any(ps.gen(:,C.ge.Pg)>ps.gen(:,C.ge.Pmax)+EPS), error('Pg is out of bounds'); end
    flow  = ps.branch(:,C.br.Pf);
    
    if nargout>5
        flows = [flows flow]; %#ok<AGROW>
    end
    % error check:
    if any(ps.gen(:,C.ge.Pg)>ps.gen(:,C.ge.Pmax)+EPS), error('Pg is out of bounds'); end

    % Step 4a. Take control actions if needed.
    if opt.sim.use_control
        % write out the status of each node in the network
        %  output file name is grid_status_{pid}.csv
        node_has_power = find_buses_with_power(ps,sub_grids);
        csvwrite(grid_status_file,[ps.bus(:,1) node_has_power]);
        
        % call python comms code
         % send the matlab pid to python code

        % check the status of the comm network
         % read the file /tmp/comm_status_{pid}.csv
         % busno,comm?
         % 1,1
         % 2,0 etc...
        if exist(comm_status_file,'file')
            data = csvread(comm_status_file);
            comm_status = data(:,2);
        else
            comm_status = true(n,1);
        end
        % Read power system data (flow) from each operating node
        is_br_readable = comm_status(F) | comm_status(T);
        measured_flow = nan(m,1);
        measured_flow(is_br_readable) = flow(is_br_readable);
        branch_st = ps.branch(:,C.br.status);
        % if there are overloads in the system, try to mitigate them
        if any(abs(measured_flow)>flow_max)
            % Use the data to choose "optimal" decisions
            ramp_dt = min(dt,dt_max); % the amount of generator ramping time to allow
            max_ramp = ramp_rate_emergency*ramp_dt;
            [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,max_ramp,opt.verbose);
            %[delta_Pd,delta_Pg] = emergency_control(ps,sub_grids,max_ramp,measured_flow,comm_status,opt.verbose);
            
            % Implement the load and generator control
            Pg_new = ps.gen(:,C.ge.P) + delta_Pg;
            ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
            % compute the new load factor
            delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
            lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
            ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));
            dt_max = 60;
            % run dcpf again
            if any(abs(delta_Pd)>EPS)
                ps = dcpf(ps,sub_grids,true);
                if any(ps.gen(:,C.ge.Pg)>ps.gen(:,C.ge.Pmax)+EPS), error('Pg is out of bounds'); end
            end
       end
    end
    % Step 5. Update relays
    [ps.relay,br_out_new,dt,n_over] = update_relays(ps,opt.verbose,dt_max);
    if opt.verbose && n_over>0
        fprintf(' There are %d overloads in the system\n',n_over);
    end
    
    % Step 6. Check for any remaining overload potential, decide if we
    % should stop the simulation
    if dt==Inf
        is_blackout = 0;
        break
    end
    % If we want to stop when the network is divided into subnetworks, do
    % this:
    if opt.sim.stop_on_sep
        if sep==BIG_SEP
            is_blackout = 1;
            if opt.verbose
                fprintf('-------------- t = %.3f ----------------\n',t);
                fprintf('----------- Major separation -----------\n');
            end
            break
        end
    else % If we wanted to stop after a certain amount of load shedding, do this:
        Pd_sum = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
        load_remaining_fraction = Pd_sum/Pd0_sum;
        if opt.verbose
            %fprintf('------------- t = %.3f ---------------\n',t);
            fprintf('-------- %.1f%% of load remains --------\n',load_remaining_fraction*100);
        end
        if load_remaining_fraction<opt.sim.stop_threshold && ~is_blackout
            is_blackout = 1;
            if opt.verbose
                fprintf('----------- Blackout occurred ----------\n');
            end
        end
        if load_remaining_fraction<opt.sim.stop_threshold
            break
        end
    end
    
    % advance/print the time
    t = t + dt;
    if opt.verbose
        fprintf('------- t = %.3f ----------\n',t);
    end
    
    % Step 7. Trip overloaded branches
    ps.branch(br_out_new,C.br.status) = 0;
    % record which branches were lost
    for i = 1:length(br_out_new)
        br = br_out_new(i);
        relay_outages = cat(1,relay_outages,[t br]);
    end
    
    % print something
    if opt.verbose && ~isempty(br_out_new)
        fprintf(' Branch %d triped on overcurrent\n',br_out_new);
    end
    
    % Return to step 3.
end

% do a final redispatch just to make sure
[Pg,ge_status,d_factor] = redispatch(ps,sub_grids,ramp_rate*dt,opt.verbose);
% error check:
if Pg>ps.gen(:,C.ge.Pmax)+EPS, error('Pg is out of bounds'); end
ps.shunt(:,C.sh.factor) = d_factor;
ps.gen(:,C.ge.status) = ge_status;
ps.gen(:,C.ge.P) = Pg;
Pd = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
MW_lost = sum(Pd0) - sum(Pd);
n_overloads = sum(ps.branch(:,C.br.Pf)>ps.branch(:,C.br.rateB));

if opt.verbose
    fprintf('-------------- t = %7.3f -----------------\n',t);
    fprintf(' Simulation complete\n');
    fprintf('  %d emergency (rateB) overloads remain\n',n_overloads);
    fprintf('  %d endogenous relay outages\n',size(relay_outages,1));
    fprintf('  %g MW load lost (%.1f%%)\n',MW_lost,MW_lost/Pd0_sum*100);
    fprintf('--------------------------------------------\n');
end

