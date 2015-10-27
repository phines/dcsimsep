function [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,ramp_limits,comm_status,opt,test_trivial)
% An emergency control function to adjust load, based on current measured
%  flows
% usage: [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,ramp_limits,comm_status,opt,test_trivial)
% inputs:
%  ps - power systems data
%  measured_flow - the amount of flow (in MW/MVA) on each transmission line
%  branch_st - the status of the transmission lines
%  opt.verbose - a binary indicating whether to print stuff about what is going
%   on.
%  
% outputs:
%  delta_Pd - changes to amount of load at each "shunt"
%  delta_Pg - changes to the generator output at each generator

%% Prep work
C = psconstants;
EPS = 1e-6;

% check the inputs
if nargin<3, error('need at least 3 inputs'); end
% number of buses in the system
n = size(ps.bus,1);
if nargin<4, ramp_limits = ps.gen(:,C.ge.Pmax); end
if nargin<5, comm_status = true(n,1); end
if nargin<6, opt.verbose = true; end
if nargin<7, test_trivial=false; end

% display something
if opt.verbose, disp('  Doing prep work for emergency control');end
% collect load data
D = ps.bus_i(ps.shunt(:,1)); % load bus locations
nd = length(D);
Pd0 = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor) ./ ps.baseMVA; % in per unit
% collect generator data
G = ps.bus_i(ps.gen(:,1));  % gen bus locations
ng = length(G);
ge_status = ps.gen(:,C.ge.status)==1;
Pg0 = ps.gen(:,C.ge.P).*ge_status./ ps.baseMVA; % in per unit
Pg = Pg0;
ramp_limits_pu = ramp_limits/ps.baseMVA;
Pg_max = Pg0; % < this could change
Pg_min = max(0,max(Pg0-ramp_limits_pu,ps.gen(:,C.ge.Pmin).*ge_status./ps.baseMVA));
if any(Pg<Pg_min-EPS) || any(Pg>Pg_max+EPS)
    error('Generation outside of Pmin/Pmax');
end
% collect transmission line data
branch_st(isnan(branch_st)) = 1; % if we don't know otherwise, assume branch is closed
br_st = (branch_st==1);
F = full(ps.bus_i(ps.branch(br_st,1)));
T = full(ps.bus_i(ps.branch(br_st,2)));
X = ps.branch(br_st,C.br.X);
inv_X = 1./X;
measured_flow_pu = measured_flow(br_st) / ps.baseMVA;
measured_flow_pu(isnan(measured_flow_pu)) = 0;
flow_max = ps.branch(br_st,C.br.rateB) / ps.baseMVA; % in per unit
m = length(F); % number of transmission lines in the system
% process the comm_status
buses = (1:n)';
comm_connected_buses = buses(comm_status);
is_G_conn = ismember(G,comm_connected_buses);
is_D_conn = ismember(D,comm_connected_buses);

%% edit this code so that it operates separately on each sub-grid
nodes = (1:n)';
links = [F,T];
[grid_no,n_sub] = find_subgraphs(nodes,links);

%% set up an index so that we can find things
ix.x.d_theta = (1:n);
ix.x.dPg = (1:ng) + n;
ix.x.dPd = (1:nd) + n + ng;
ix.x.f_over = (1:m) + n + ng + nd;
ix.nx = n + ng + nd + m;

%% set up limits and costs for the variables
x_min = zeros(ix.nx,1)-Inf;
x_max = zeros(ix.nx,1)+Inf;
% for dPg
x_min(ix.x.dPg) = min(Pg_min - Pg0,0).*is_G_conn;
x_max(ix.x.dPg) = (Pg_max - Pg0).*is_G_conn;
% for dPd
x_min(ix.x.dPd) = min(-Pd0,0).*is_D_conn;
x_max(ix.x.dPd) = 0;
% for f_over
x_min(ix.x.f_over) = 0;
x_max(ix.x.f_over) = Inf;
% constrain one reference bus for each island
% choose a bus that is already closest to zero
theta = ps.bus(:,C.bu.Vang) * pi / 180;
ref_bus_num = zeros(n_sub,1);
for grid_i = 1:n_sub
    bus_subset = find(grid_no == grid_i);
    theta_sub = theta(bus_subset);
    % find the bus in this island with the smallest absolute angle
    [~,ref_tmp] = min(abs(theta_sub));
    ref_ix = bus_subset(ref_tmp);
    % Add this bus to the reference list
    ref_bus_num(grid_i) = ref_ix;
end
% set this bus to be the local reference
x_min(ix.x.d_theta(ref_bus_num)) = 0;
x_max(ix.x.d_theta(ref_bus_num)) = 0;

cost = zeros(ix.nx,1);
cost(ix.x.dPd) = -opt.sim.cost.load;
cost(ix.x.f_over) = opt.sim.cost.overload;

%% set up equality constraints
% DC power flow constraint
A_pf = sparse(n,ix.nx); 
% generator injection constraint:
A_pf = A_pf + sparse(G,ix.x.dPg,1,n,ix.nx);
% load injection constraint:
A_pf = A_pf + sparse(D,ix.x.dPd,-1,n,ix.nx);
% DC power flow matrix constraint
A_pf = A_pf + sparse(F,T,+inv_X,n,ix.nx) + ...
              sparse(T,F,+inv_X,n,ix.nx) + ...
              sparse(T,T,-inv_X,n,ix.nx) + ...
              sparse(F,F,-inv_X,n,ix.nx);
     
b_pf = zeros(n,1);

%% set up the inequality constraints
% These are the two flow constraints, which are:
% (d_theta_f-d_theta_t)*1/Xft - f_over <= -f_meas + f_max   (1)
% (d_theta_f-d_theta_t)*1/Xft + f_over >= -f_meas - f_max   (2)
A_flow = sparse((1:m),F,+inv_X,m,ix.nx) + ... % Implements 1/Xft for the F side of the lines
         sparse((1:m),T,-inv_X,m,ix.nx);     % Implements 1/Xft for the T side of the lines
% add the stuff for Eq (1)
A_flow_1 = +A_flow + sparse(1:m,ix.x.f_over,-1,m,ix.nx);
% add the stuff for Eq (2)
A_flow_2 = +A_flow + sparse(1:m,ix.x.f_over,+1,m,ix.nx);
% implement the right side of these equations
b_flow_1 = -measured_flow_pu + flow_max;
b_flow_2 = -measured_flow_pu - flow_max;
A_ineq = [A_flow_1;-A_flow_2];
b_ineq = [b_flow_1;-b_flow_2];

%% if requested test the trivial solution
if test_trivial
    x_test = zeros(ix.nx,1);
    x_test(ix.x.d_theta) = 0;
    x_test(ix.x.dPd) = 0;
    x_test(ix.x.dPg) = 0;
    x_test(ix.x.f_over) = abs(measured_flow_pu);
    disp('X bounds');
    violation = ~(x_min-EPS<x_test & x_test<x_max+EPS);
    %violation = ~(x_test<x_max+EPS);
    disp([x_min x_test x_max violation])
    if any(violation)
        error('Trivial solution doesn''t work');
    end
    disp('Equality cons');
    violation = (abs(A_pf*x_test - b_pf)>EPS);
    disp([A_pf*x_test - b_pf violation]);
    if any(violation)
        error('Trivial solution doesn''t work');
    end
    disp('Inequality cons');
    Ax = A_ineq*x_test;
    violation = ~(Ax<=(b_ineq+EPS));
    disp([Ax b_ineq violation]);
    if any(violation)
        error('Trivial solution doesn''t work');
    end
end

%% run the optimization
if opt.verbose
    disp('  Solving the emergency control problem');
end
switch opt.optimizer
    case 'gurobi'
        [x_star, fval, exitflag] = gurobi_lp(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);
    case 'cplex'
        [x_star,fval,exitflag,output] = cplexlp(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);
    case 'linprog'
        [x_star,fval,exitflag,output] = linprog(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);
    
    %{
        case 'mexosi'
        error('mexosi is broken');
        A = [A_pf;A_flow_1;A_flow_2];
        b_max = [b_pf;b_flow_1;b_flow_2];
        b_min = -Inf(size(b_max));
        [x_star,~,exitflag] = osi(cost,x_min,x_max,A,b_min,b_max);
        %[x_star,~,exitflag,~] = linprog(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);
        %}
    case 'cvx'
        B = sparse(F,T,-inv_X,n,n) + ...
            sparse(T,F,-inv_X,n,n) + ...
            sparse(T,T,+inv_X,n,n) + ...
            sparse(F,F,+inv_X,n,n);
        A_Pg = zeros(n,ng); IDX = sub2ind(size(A_Pg),G',1:length(G)); A_Pg(IDX) = 1;
        A_Pd = zeros(n,nd); IDX = sub2ind(size(A_Pd),D',1:length(D)); A_Pd(IDX) = 1;
        warning('off','MATLAB:nearlySingularMatrix'); % suppress this warning
        if opt.verbose 
            cvx_begin
        else
            cvx_begin quiet
        end
        variables delta_Pg_pu(ng) delta_Pd_pu(nd) f_over(m) delta_theta(n)
        minimize(-sum(delta_Pd_pu) + f_over_cost*ones(1,m)*f_over)
        subject to:
            -flow_max - f_over <= measured_flow_pu + diag(inv_X)*(delta_theta(F)-delta_theta(T)) <= flow_max + f_over;
            B*delta_theta == A_Pg*delta_Pg_pu - A_Pd*delta_Pd_pu;
            min(Pg_min - Pg0,0).*is_G_conn <= delta_Pg_pu <= (Pg_max - Pg0).*is_G_conn;
            min(-Pd0,0).*is_D_conn <= delta_Pd_pu <= 0;
            f_over >= 0;
            delta_theta(ref_bus_num) == 0;
        cvx_end
        if strcmp(cvx_status,'Solved')
            exitflag = 1000;
        else
            exitflag = 0;
        end
        fval = cvx_optval;
        output = cvx_status;
        warning('on','MATLAB:nearlySingularMatrix'); % turn the warning back on
end
%% check and process the solution
if exitflag==1
    delta_Pd_pu = x_star(ix.x.dPd);
    delta_Pg_pu = x_star(ix.x.dPg);
end
if exitflag==1 || exitflag==1000
    delta_Pg = delta_Pg_pu * ps.baseMVA;
    delta_Pd = delta_Pd_pu * ps.baseMVA;
    if opt.verbose
        fprintf('  Solved the emergency control problem. fval = %g.\n',fval);
    end
else
    error('  Optimization failed');
    delta_Pg_pu = zeros(ng,1);
    delta_Pd_pu = zeros(nd,1);
    delta_Pg = zeros(ng,1);
    delta_Pd = zeros(nd,1);
end

%% Debug stuff
%double check that the result is within bounds
Pg = Pg0 + delta_Pg_pu;
if any( Pg+EPS < Pg_min | Pg > Pg_max+EPS )
    error('Control results are outside of generator limits');
    %p = find(Pg+EPS < Pg_min | Pg > Pg_max+EPS);
    %[ x_min(ix.x.dPg(p)) Pg_min(p) delta_Pg(p) Pg0(p) Pg(p) Pg_max(p)]
end
Pd = Pd0 + delta_Pd_pu;
if any( Pd0>0 & ( Pd+EPS < 0 | Pd > Pd0+EPS ) )
    error('Control results are outside of load limits');
    %p = find( Pd0>0 & ( Pd < 0 | Pd > Pd0 ) )
end
%flow_pu = measured_flow_pu + A_flow * x_star;

