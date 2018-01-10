
load texas10K
ps = updateps(ps);
ps = dcpf(ps)
C = psconstants;
n = size(ps.bus,1);
n_br = size(ps.branch,1);
opt = psoptions;

br_outages = find(rand(n_br,1)<0.005);
length(br_outages)

opt.verbose = 1;
opt.sim.stop_threshold = 0.95

is_bo = dcsimsep(ps,br_outages,[],opt);

%printps(ps)
