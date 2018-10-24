function opt = psoptions(opt)
% some options for the power system simulator files

if nargin==0
    opt = struct;
end

C = psconstants;

%% General options
opt.verbose = 1;
%opt.seecascade = 0;
opt.debug = true;

opt.optimizer = 'gurobi'; % or cplex, mexosi

%% power flow options
opt.pf.tolerance = 1e-8; % convergence tolerance
opt.pf.max_iters = 20; % max power flow iterations
opt.pf.CalcIslands = 1; % iteratively calculates each island in runPowerFlow
opt.pf.CascadingPowerFlow = 0;
opt.pf.flat_start = 0;
opt.pf.load_shed_rate = 0.25; % the rate at which under frequency load shedding is done in CascadingPowerFlow mode
opt.pf.linesearch = 'exact';
opt.pf.update = true;
opt.pf.check_Pg = true;

%% optimal power flow options
opt.opf.generator_commitment = 0; % switch generators on/off using MIP
opt.opf.branch_switching = 0;     % switch branches on/off using MIP
opt.opf.rate = C.br.rateA;
opt.opf.contg_rate = C.br.rateB;

%% time-domain simulation options options
opt.sim.ramp_frac = 0.05; % fraction of generator allowed to ramp between generations
opt.sim.writelog  = false;
% opt.sim.dt_default = 10; % default (max) time step size  ----->  use opt.sim.dt instead of this 
opt.sim.dt = 60; % maximum time step between cascade simulation iterations
opt.sim.t_max = 15*60; % maximum simulation time in dcsimsep
opt.sim.draw = true;
opt.sim.overload_time_limit = 10*60; % number of seconds that the branch can sit at its rateC level (PSS/E manual)
opt.sim.stop_on_sep = false; % When true, the simulator declares a blackout
                             % based on the Giant Component size. Otherwise
                             % this depends on the blackout size in MW
opt.sim.stop_threshold = 0.95; % the fraction of nodes in the giant component, or load still connected, at which to stop the simulation
opt.sim.fast_ramp_mins = 1; % the minimum minutes of ramping that generators are allowed to do without load shedding
% opt.sim.use_control = false; -----> use opt.sim.control_method instead
opt.sim.control_method = 'none';
opt.sim.use_comm_model = false;
% opt.sim.dt_max_default = 60; % maximum amount of time between dcsimsep iterations   ----->  use opt.sim.dt instead of this 
opt.sim.simple_rebalance = false; % Simple method used by Zussman's model
opt.sim.relay_trip_time = 15; % time to trip an overcurrent relay when the line is at rateC
opt.sim.cost.load = 1;
opt.sim.cost.overload = 1000;
% distributed control options
opt.sim.nHopLoc = 2;
opt.sim.nHopExt = 10;
% MPC options
opt.sim.Np = 3;
opt.sim.use_mpc = false;


%% drawing stuff
opt.draw.width = 0.1;
opt.draw.bus_nos = true;
opt.draw.simple = false;
opt.draw.fontsize = 14;


