% test simulate_dc
clear all;

C = psconstants; % tells me where to find my data
opt = psoptions;
opt.verbose = true;
opt.sim.stop_threshold = 0.90; % the fraction of nodes, at which to declare a major separation

%% a simple test
%{
disp('Running 6-bus test');
ps = updateps(case6ww_ps);
ps.gen(:,C.ge.ramp_rate_up) = ps.gen(:,C.ge.Pmax)/10; %10% per minute
ps = dcpf(ps);
m = size(ps.branch,1);
% choose some outages
br_outages = randi(m,2,1);
% run the simulator
[is_blackout,relay_outages,MW_lost] = dcsimsep(ps,br_outages,[],opt);
if is_blackout
    disp('Blackout');
else
    disp('Not blackout');
end
MW_lost
%}

%% a more complicated test
disp('Running 2383 bus test');

% input parameters:
p_outage = 0.005;
%casename = 'case300_001_ps');
disp('loading the data');
tic
%casename = 'case2383_mod_ps';
%ps = updateps(feval(casename));
load case2383_mod_ps;
m = size(ps.branch,1);
toc
opt.verbose = 1;

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


