clear all
%% create the case
C = psconstants;
opt = psoptions;
ps = case6ww_ps;
ps = updateps(ps);
ps = dcpf(ps);
printps(ps);
ramp_limits = [100;100;100];
verbose = 1;

%% collect some data
n = size(ps.bus,1);
m = size(ps.branch,1);
ng = size(ps.gen,1);
nd = size(ps.shunt,1);
Pg_min = zeros(ng,1);
Pg_max = ps.gen(:,C.ge.Pmax);

%% test split the network
% split the network in two
ps.branch([1 5 6 8 11],C.br.status) = 0;
% run redispatch
[sep,sub_grids,n_sub] = check_separation(ps,opt.sim.stop_threshold,opt.verbose);
ps = redispatch(ps,sub_grids,ramp_limits,verbose);
ps = dcpf(ps);
printps(ps);

%% measure the network state
measured_flow = ps.branch(:,C.br.Pf);
branch_st = ps.branch(:,C.br.status);
printps(ps);
figure(1);
drawps(ps,opt);

%% calculate the control
verbose = true;
[delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,verbose);

%% implement the control
% Implement the load and generator control
Pg_new = ps.gen(:,C.ge.P) + delta_Pg;
ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
% compute the new load factor
delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));

%% run power flow and draw
ps = dcpf(ps);
figure(2);
drawps(ps,opt);
measured_flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
[measured_flow flow_max]

%% repeat the test for the 30 bus case
C = psconstants;
opt = psoptions;
ps = case30_ps;
ps = updateps(ps);
ps = dcpf(ps);
printps(ps);
verbose = 1;
% collect some data
n = size(ps.bus,1);
m = size(ps.branch,1);
ng = size(ps.gen,1);
nd = size(ps.shunt,1);
Pg_min = zeros(ng,1);
Pg_max = ps.gen(:,C.ge.Pmax);
ramp_limits = Pg_max;
% edit the line limits
ps.branch(:,C.br.rates) = ps.branch(:,C.br.rates)*.5;
% remove lines
% split the network in three
ps.branch([33 36],C.br.status) = 0;
ps.branch([12 14 15],C.br.status) = 0;
% run redispatch
[sep,sub_grids,n_sub] = check_separation(ps,opt.sim.stop_threshold,opt.verbose);
ps = redispatch(ps,sub_grids,ramp_limits,verbose);
ps = dcpf(ps);
% measure the network state
measured_flow = ps.branch(:,C.br.Pf);
branch_st = ps.branch(:,C.br.status);
printps(ps);
figure(3);
drawps(ps,opt);

% calculate the control
verbose = true;
[delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,verbose,1);
% implement the control
% Implement the load and generator control
Pg_new = ps.gen(:,C.ge.P) + delta_Pg;
ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
% compute the new load factor
delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));
% run power flow and draw
ps = dcpf(ps);
figure(4);
drawps(ps,opt);
measured_flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
[measured_flow flow_max]

