% test simulate_dc
clear all;

C = psconstants; % tells me where to find my data
opt = psoptions;
opt.verbose = false;
opt.sim.stop_threshold = 0.90; % the fraction of nodes, at which to declare a major separation

%% a more complicated test
disp('Testing each br outage for the 2383 bus test');

% load the data
disp('loading the data');
load case2383_mod_ps;
ps_org = ps;

% take out each branch
m = size(ps.branch,1);
for i = 1:m
    ps = ps_org;
    br_out = i;
    [is_blackout(i),relay_outages,MW_lost(i)] = dcsimsep(ps,br_out,[],opt);
    if MW_lost(i)>50
        fprintf('Br: %d, BO: %d, MW lost: %g\n',i,is_blackout(i),MW_lost(i));
    end
end

return

%% run a set of simulations
opt.verbose = 0;
% choose some branch outages
bus_outages = [];
% run the simulator
for i=1:1000
    br_outages = find(rand(m,1)<p_outage)';
    is_blackout(i) = dcsimsep(ps,br_outages,[],opt);
    if is_blackout(i)
        fprintf('%5d: %2d outages. Blackout\n',i,length(br_outages));
    else
        fprintf('%5d: %2d outages. Not blackout\n',i,length(br_outages));
    end
end


