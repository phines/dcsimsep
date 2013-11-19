
%% create the case
C = psconstants;
ps = case6ww_ps;
ps = updateps(ps);
ps = dcpf(ps);
measured_flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
ps.branch(2,C.br.rates) = 40;
printps(ps);
opt = psoptions;
figure(1);
drawps(ps,opt);
pause

%% calculate the control
verbose = true;
[delta_Pd,delta_Pg] = emergency_control_simple(ps,measured_flow,verbose);

%% implement the control
ps.shunt(:,C.sh.P) = ps.shunt(:,C.sh.P) - delta_Pd;
ps.gen(:,C.ge.P) = ps.shunt(:,C.ge.P) - delta_Pg;
ps = dcpf(ps);
figure(2);
drawps(ps,opt);
pause
%[measured_flow flow_max]
