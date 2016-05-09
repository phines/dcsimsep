% test simulate_dc
clear all;

% constants
C = psconstants; % tells me where to find my data
opt = psoptions;
opt.verbose = true;
opt.sim.stop_threshold = 0.90; % the fraction of nodes, at which to declare a major separation

%% load the data
load case2383_mod_ps;
ps = case6_ps;
ps = updateps(ps);
ps = rebalance(ps);
ps = dcpf(ps);

%% collect some data about the system
m = size(ps.branch,1);
flows = ps.branch(:,C.br.Pf) / ps.baseMVA;
flow_max = ps.branch(:,C.br.rateA) / ps.baseMVA;

% build the distribution factors
PTDF = makePTDF(ps.baseMVA,ps.bus,ps.branch);
LODF = makeLODF(ps.branch,PTDF);
% get rid of the nans in the lodf
LODF(isnan(LODF)) = 0;

% calculate the matrix of post-contingency flows
flows_prime = ones(m,1)*flows' + LODF * diag(flows);
is_over_limit = flows_prime > ones(m,1)*flow_max';
flows_pp    = flows_prime .* is_over_limit;
vulnerability_metric = sum(flows_pp); % is this useful???

return

return


% choose some branch outages
%br_outages = [17 24 51 143 365];
%br_outages = find(rand(m,1)<p_outage)';
br_outages = [152];
bus_outages = [];

% run the simulator
disp('running the simulation');
[is_blackout,relay_outages,MW_lost] = dcsimsep(ps,br_outages,bus_outages,opt);

if is_blackout
    disp('Blackout');
else
    disp('Not blackout');
end

% draw the cascade
if is_blackout
    draw_cascade(ps,br_outages,bus_outages,relay_outages);
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


