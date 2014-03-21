clear all;

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation
opt.sim.fast_ramp_mins = 1;

%% Prepare and run the simulation for the 6 bus case
% load the case data
ps = case6ww_ps;
ps = updateps(ps);
ps = redispatch(ps);
ps = dcpf(ps);
%printps(ps);
verbose = 1;
EPS = 1e-6;
% collect some data about the case
n = size(ps.bus,1);
m = size(ps.branch,1);
ng = size(ps.gen,1);
nd = size(ps.shunt,1);
Pg_min = zeros(ng,1);
Pg_max = ps.gen(:,C.ge.Pmax);
flow_max = ps.branch(:,C.br.rateB);

% choose some branch outages
br_outages = [1 5 6 8 11];
%br_outages = [2 4];
% run the simulator
[is_blackout,outages,MW_lost] = dcsimsep(ps,br_outages,[],opt)

%% try again with emergency control
opt.sim.use_control = 1;
% run the simulator
[is_blackout,outages,MW_lost] = dcsimsep(ps,br_outages,[],opt)


return

