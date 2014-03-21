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
ps = redispatch(ps);
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
flows_prime = flows*ones(m,1)' + LODF * diag(flows);
% check to see if these exceed the thresholds
is_over_limit = abs(flows_prime) > flow_max*ones(m,1)';
% weed out the non threshold crossings
flows_pp    = flows_prime .* is_over_limit;
% sum up the flows on the lines after exceeding the thresholds
vulnerability_metric = sum(abs(flows_pp),1) % is this useful???

return
