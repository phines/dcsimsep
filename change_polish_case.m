% test simulate_dc
clear all;

C = psconstants; % tells me where to find my data
opt = psoptions;
opt.verbose = true;
opt.sim.stop_threshold = 0.95; % the fraction of nodes, at which to declare a major separation

%% a more complicated test
disp('Changing the polish test case');

casename = 'case2383_mod_ps';
ps = updateps(feval(casename));

X = ps.branch(:,C.br.X);
rateA = 1./X * ps.baseMVA / 100;
rateB = rateA * 1.2;
rateC = rateA * 1.5;

ps.branch(:,C.br.rateA) = rateA;
ps.branch(:,C.br.rateB) = rateB;
ps.branch(:,C.br.rateC) = rateC;


%% run a simulation

% input parameters:
p_outage = 0.005;
m = size(ps.branch,1);
opt.verbose = 1;

% choose some branch outages
%br_outages = [17 24 51 143 365];
br_outages = find(rand(m,1)<p_outage)';
%br_outages = [96];
bus_outages = [];

% run the simulator
disp('running the simulation');
[is_blackout,relay_outages] = dcsimsep(ps,br_outages,[],opt);

if is_blackout
    disp('Blackout');
else
    disp('Not blackout');
end

% draw the cascade
draw_cascade(ps,br_outages,bus_outages,relay_outages);
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


