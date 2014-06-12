clear all
%% Global constants
C = psconstants;
opt = psoptions;
opt.verbose = 1;
verbose = opt.verbose;
opt.optimizer = 'cplex';
EPS = 1e-6;

%% create the case
ps = case6ww_ps;
ps = updateps(ps);
ps = dcpf(ps);
printps(ps);

%% collect some data
n = size(ps.bus,1);
m = size(ps.branch,1);
ng = size(ps.gen,1);
nd = size(ps.shunt,1);
Pg_min = zeros(ng,1);
Pg_max = ps.gen(:,C.ge.Pmax);
flow_max = ps.branch(:,C.br.rateB);
comm_status = true(n,1);

%% test split the network
% split the network in two
ps.branch([1 5 6 8 11],C.br.status) = 0;
% run redispatch
[sep,sub_grids,n_sub] = check_separation(ps,opt.sim.stop_threshold,opt.verbose);
ramp_limits = Pg_max;
ps = redispatch(ps,sub_grids,ramp_limits,verbose);
verbose = true;
ramp_limits = Pg_max*.1;

%% iteratively calculate and implement the control
while 1
    % run power flow
    ps = dcpf(ps);
    % collect the state data
    branch_st = ps.branch(:,C.br.status);
    measured_flow = ps.branch(:,C.br.Pf);
    if all( measured_flow <= (flow_max+EPS) )
        break;
    end
    figure(1);
    drawps(ps,opt);
    title('Before');
    pause
    % optimize
    [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,ramp_limits,comm_status,opt)
    % implement the control
    % Implement the load and generator control
    Pg_new = ps.gen(:,C.ge.P) + delta_Pg;
    ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
    % compute the new load factor
    delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
    lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
    ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));
end

%% run power flow and draw
ps = dcpf(ps);
figure(2);
drawps(ps,opt);
measured_flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
[measured_flow flow_max]
title('After');
pause

%% repeat the test for the 30 bus case
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
flow_max = ps.branch(:,C.br.rateB);
% remove lines
% split the network in three
ps.branch([33 36],C.br.status) = 0;
ps.branch([12 14 15],C.br.status) = 0;
% run redispatch
[sep,sub_grids,n_sub] = check_separation(ps,opt.sim.stop_threshold,opt.verbose);
ps = redispatch(ps,sub_grids,ramp_limits,verbose);
% iteratively calculate and implement the control
ramp_limits = Pg_max * 0.1;
comm_status = true(n,1);

while 1
    % run power flow
    ps = dcpf(ps);
    % collect the state data
    branch_st = ps.branch(:,C.br.status);
    measured_flow = ps.branch(:,C.br.Pf);
    if all( measured_flow <= (flow_max+EPS) )
        break;
    end
    figure(3);
    drawps(ps,opt);
    title('Before');
    pause
    % optimize
    [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,ramp_limits,comm_status,opt)
    % implement the control
    % Implement the load and generator control
    Pg_new = ps.gen(:,C.ge.P) + delta_Pg;
    ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
    % compute the new load factor
    delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
    lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
    ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));
end

%% run power flow and draw
ps = dcpf(ps);
figure(4);
drawps(ps,opt);
measured_flow = ps.branch(:,C.br.Pf);
title('After');

%% Now try with the Polish case
disp('Checking the Polish case');
ps = case2383_mod_ps;
ps = updateps(ps);
verbose = 0;
ps = redispatch(ps,[],[],verbose);
ps = dcpf(ps);
ps = redispatch(ps,verbose);
% collect some data
n = size(ps.bus,1);
m = size(ps.branch,1);
ng = size(ps.gen,1);
nd = size(ps.shunt,1);
Pg_min = zeros(ng,1);
Pg_max = ps.gen(:,C.ge.Pmax);
ramp_limits = Pg_max;
% edit the line limits if wanted
%ps.branch(:,C.br.rates) = ps.branch(:,C.br.rates)*.5;
flow_max = ps.branch(:,C.br.rateB);
% remove lines
% choose a contingency
ps.branch([11 38],C.br.status) = 0;
% run redispatch
[sep,sub_grids,n_sub] = check_separation(ps,opt.sim.stop_threshold,opt.verbose);
ps = redispatch(ps,sub_grids,ramp_limits,verbose);
% iteratively calculate and implement the control
ramp_limits = Pg_max * 0.5;
while 1
    % run power flow
    ps = dcpf(ps);
    % collect the state data
    branch_st = ps.branch(:,C.br.status);
    measured_flow = ps.branch(:,C.br.Pf);
    if all( abs(measured_flow) <= (flow_max+EPS) )
        break;
    end
    % optimize
    comm_status = true(n,1);
    [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,ramp_limits,comm_status,opt);
    % Implement the load and generator control
    Pg_new = ps.gen(:,C.ge.P) + delta_Pg;
    ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
    % compute the new load factor
    delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
    lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
    ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));
end

measured_flow = ps.branch(:,C.br.Pf);
n_violations = sum(abs(measured_flow) > flow_max+EPS)

