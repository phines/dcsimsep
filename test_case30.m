clear all;

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.0; % the threshold at which to declare a blackout
opt.verbose = 1;
opt.sim.use_control = 1;

%% Prepare and run the simulation for the Polish grid
disp('loading the data');
ps = case30_ps;
ps = updateps(ps);
ps = dcpf(ps);
m = size(ps.branch,1);
ps.gen

% choose a set of component (transmission line) outages. 
% For the 30-bus case this can be any set of integers from 1 to 41
br_outages = [23 6 7];
bus_outages = [];

% run the simulator
disp('running the simulation');
[is_blackout,relay_outages,MW_lost] = dcsimsep(ps,br_outages,bus_outages,opt);

%
if is_blackout
    disp('Blackout');
else
    disp('Not blackout');
end

%% draw the cascade
if is_blackout
    draw_cascade(ps,br_outages,bus_outages,relay_outages);
end

return

