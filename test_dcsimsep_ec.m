clear all;

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation
opt.sim.fast_ramp_mins = 1;
opt.sim.use_control = 1;

%% Prepare and run the simulation for the Polish grid
%ps = case300_001_ps;
fprintf('----------------------------------------------------------\n');
disp('loading the data');
tic
load case2383_mod_ps;
toc
fprintf('----------------------------------------------------------\n');
tic
ps = updateps(ps);
ps = dcpf(ps);
toc
fprintf('----------------------------------------------------------\n');
m = size(ps.branch,1);
pre_contingency_flows = ps.branch(:,C.br.Pf);
phase_angles_degrees = ps.bus(:,C.bu.Vang);

% choose some branch outages
load BOpairs
n_iters = 100;
tic
for i = 1:n_iters
    % outage
    br_outages = BOpairs(i,:);
    % run the simulator
    fprintf('Running simulation %d of %d\n',i,n_iters);
    %[is_blackout(i),relay_outages,MW_lost(i),p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);
    is_blackout = dcsimsep(ps,br_outages,[],opt);
end
toc

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

