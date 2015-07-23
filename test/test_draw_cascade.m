% test simulate_dc
clear all;
% constants
C = psconstants; % tells me where to find my data
opt = psoptions;
opt.verbose = true;
opt.sim.stop_threshold = 0.90; % the fraction of nodes, at which to declare a major separation

% input data
sim_set = [1:10];

%% run the simulation
% load the data
disp('loading the data');
load case2383_mod_ps;
load BOpairs;
ps_org = ps;
pairs = BOpairs;
% take out each branch
m = size(ps.branch,1);
for i = sim_set
    ps = ps_org;
    brs_out = pairs(i,:);
    [is_blackout,relay_outages,MW_lost,p_out,busessep] = dcsimsep(ps,brs_out,[],opt);
    draw_cascade(ps,brs_out,[],relay_outages);
    pause
end

return
