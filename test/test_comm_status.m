clear all;
clc

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation
opt.sim.fast_ramp_mins = 1;
opt.sim.use_control = true;
opt.sim.use_comm_model = true;
opt.comm.two_way = true;

%% Prepare the simulation for the Polish grid
fprintf('----------------------------------------------------------\n');
disp('loading the data');
tic
if exist('case2383_mod_ps.mat','file')
    load case2383_mod_ps;
else
    ps = case2383_mod_ps;
end
toc
fprintf('----------------------------------------------------------\n');
tic
ps = updateps(ps);
ps = redispatch(ps);
ps = dcpf(ps);
toc
fprintf('----------------------------------------------------------\n');
m = size(ps.branch,1);
pre_contingency_flows = ps.branch(:,C.br.Pf);
phase_angles_degrees = ps.bus(:,C.bu.Vang);
% Assign each bus to exactly one load (for station service)
ps.bus(:,C.bu.power_from_sh) = assign_loads_to_buses(ps);

%% write out a comm status file.
pid = feature('getpid');
cmd = sprintf('cp ./comm_status_test.csv /tmp/comm_status_%d.csv',pid);
system(cmd);

%% run a single case
load BOpairs
opt.verbose = true;
i = 247;
br_outages = BOpairs(i,:);
[~,relay_outages,MW_lost_1(i),p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);

return
%% Run several cases
opt.verbose = false;

n_iters = 50;
% try again with control
opt.sim.use_control = false;
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost_1(i),p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt); %#ok<SAGROW>
    fprintf(' Result: %.2f MW of load shedding\n',MW_lost_1(i));
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

% try again with control
opt.sim.use_control = true;
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages_2,MW_lost_2(i),p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt); %#ok<SAGROW>
    fprintf(' Result: %.2f MW of load shedding\n',MW_lost_2(i));
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

return
%% make a picture
figure(1);
Pd0 = sum(ps.shunt(:,C.sh.P));
bo_sizes = MW_lost/Pd0;
bar(1:n_iters,bo_sizes,1,'EdgeColor','none');

fprintf('%d of %d simulations resulted in 10%% or more load shedding\n',sum(bo_sizes>0.1),n_iters);
return
%% draw the cascade
if is_blackout
    draw_cascade(ps,br_outages,bus_outages,relay_outages);
end

return

