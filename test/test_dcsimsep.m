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
opt.sim.simple_redispatch = false;

%% Prepare and run the simulation for the Polish grid
%ps = case300_001_ps;
fprintf('----------------------------------------------------------\n');
disp('loading the data');
tic
% if exist('case2383_mod_ps.mat','file')
%     load case2383_mod_ps;
% else
%     ps = case2383_mod_ps;
% end
%load ps_polish_100;
%ps = ps_polish_100;
load case2383_mod_ps
toc
fprintf('----------------------------------------------------------\n');
tic
ps = updateps(ps);
ps = rebalance(ps);
ps = dcpf(ps);
toc
fprintf('----------------------------------------------------------\n');
m = size(ps.branch,1);
pre_contingency_flows = ps.branch(:,C.br.Pf);
phase_angles_degrees = ps.bus(:,C.bu.Vang);
Pd_total = sum(ps.shunt(:,C.sh.P));
% Set lower gen limits to zero
ps.gen(:,C.ge.Pmin) = 0;

%% test one case of interest
br_outages = [23 168 169];
[~,relay_outages,MW_lost] = dcsimsep(ps,br_outages,[],opt);
total_lost = MW_lost.rebalance + MW_lost.control

%% test writeps
writeps(ps,'temp_ps');
ps = updateps(temp_ps);
ps = dcpf(ps);
printps(ps);

return

%% Run several large blackout cases
opt.verbose = false;
load BOpairs
n_iters = 10;

disp('Testing DCSIMSEP without control.');
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost] = dcsimsep(ps,br_outages,[],opt);
    total_lost = MW_lost.rebalance + MW_lost.control;
    fprintf(' Result: %.2f MW of load shedding\n',total_lost);
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

return 

%% try again with control
disp('Testing DCSIMSEP with control.');
disp('Note that this only works if you have cplex installed');
opt.sim.use_control = true;
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages_2,MW_lost_2(i),p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    fprintf(' Result: %.2f MW of load shedding\n',MW_lost_2(i));
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end

return

