clear all;
%addpath('../dcsimsep/');

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
% opt.sim.control_method = 'none';
opt.sim.control_method = 'distributed_control';
opt.sim.nHopExt = 1;
opt.sim.nHopLoc = 1;
% opt.sim.control_method = 'emergency_control';
opt.verbose = true;
opt.sim.t_max = 30*60;

%% Prepare and run the simulation for the Polish grid
%ps = case300_001_ps;
fprintf('----------------------------------------------------------\n');
disp('loading the data');
% tic
% if exist('case2383_mod_ps.mat','file')
%     load case2383_mod_ps;
% else
%     ps = case2383_mod_ps;
% end
% toc
load('ps_polish_all','ps_polish_100');
ps = ps_polish_100;
%ps = case30_mod_ps;

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
% set negative loads to zero (in Polish system)
ps.shunt(ps.shunt(:,C.sh.P)<0,C.sh.P) = 0;
ps = rebalance(ps,[],[],opt);

%% Run one extreme case
%load crashedINsimulateDC br_outages;
outage_number = 12;
load ../data/BOpairs2;
br_outages = BOpairs(outage_number,:);
% opt.control_method = 'emergency_control';
[is_bo,~,MW_lost] = dcsimsep(ps,br_outages,[],opt);
return

%% Run some extreme cases
opt.verbose = false;
n_iters = 100;
for i = 1:n_iters
    br_outages = choose_k(1:m,80);
    fprintf('Running extreme case %d of %d. ',i,n_iters);
    [is_bo,~,MW_lost] = dcsimsep(ps,br_outages,[],opt);
    fprintf(' Result: %.2f MW (%.2f%%) of load shedding\n',MW_lost,MW_lost/Pd_total*100);
end

%% Run several cases
opt.verbose = false;
n_iters = 20;

disp('Testing DCSIMSEP without control.');
opt.sim.control_method = 'none';
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    fprintf([' load shedding result: %.2f MW in rebalance and %.2f MW in ',...
        'control (%.2f%% total)\n'], MW_lost.rebalance, MW_lost.control,...
        (MW_lost.rebalance + MW_lost.control)/Pd_total*100);
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

% try again with control
fprintf('\nTesting DCSIMSEP with emergency control.\n');
opt.sim.control_method = 'emergency_control';
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    fprintf([' load shedding result: %.2f MW in rebalance and %.2f MW in ',...
        'control (%.2f%% total)\n'], MW_lost.rebalance, MW_lost.control,...
        (MW_lost.rebalance + MW_lost.control)/Pd_total*100);
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

% try again with control
fprintf('\nTesting DCSIMSEP with distributed control.\n');
opt.sim.control_method = 'distributed_control';
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    fprintf([' load shedding result: %.2f MW in rebalance and %.2f MW in ',...
        'control (%.2f%% total)\n'], MW_lost.rebalance, MW_lost.control,...
        (MW_lost.rebalance + MW_lost.control)/Pd_total*100);
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

% try again with control
fprintf('\nTesting DCSIMSEP with distributed MPC (N=1, just to compare with the above).\n');
opt.sim.control_method = 'distributed_control';
opt.sim.use_mpc = 1;
opt.sim.Np = 1;
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    fprintf([' load shedding result: %.2f MW in rebalance and %.2f MW in ',...
        'control (%.2f%% total)\n'], MW_lost.rebalance, MW_lost.control,...
        (MW_lost.rebalance + MW_lost.control)/Pd_total*100);
    %is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

% try again with control
fprintf('\nTesting DCSIMSEP with distributed MPC (N=5).\n');
opt.sim.control_method = 'distributed_control';
opt.sim.use_mpc = 1;
opt.sim.Np = 5;
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d. ',i,n_iters);
    [~,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    fprintf([' load shedding result: %.2f MW in rebalance and %.2f MW in ',...
        'control (%.2f%% total)\n'], MW_lost.rebalance, MW_lost.control,...
        (MW_lost.rebalance + MW_lost.control)/Pd_total*100);
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

