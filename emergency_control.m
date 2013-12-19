function [delta_Pd,delta_Pg] = emergency_control(ps,measured_flow,branch_st,ramp_limits,comm_status,verbose,test_trivial)
% An emergency control function to adjust load, based on current measured
%  flows
% usage: [delta_df,delta_Pg] = emergency_opf(ps,measured_flow,verbose)
% inputs:
%  ps - power systems data
%  measured_flow - the amount of flow (in MW/MVA) on each transmission line
%  branch_st - the status of the transmission lines
%  verbose - a binary indicating whether to print stuff about what is going
%   on.
%  
% outputs:
%  delta_Pd - changes to amount of load at each "shunt"
%  delta_Pg - changes to the generator output at each generator

%% Prep work
C = psconstants;
EPS = 1e-6;
f_over_cost = 1000;

% check the inputs
if nargin<3, error('need at least 3 inputs'); end
% number of buses in the system
n = size(ps.bus,1);
if nargin<4, ramp_limits = ps.gen(:,C.ge.Pmax); end
if nargin<5, comm_status = true(n,1); end
if nargin<6, verbose = true; end
if nargin<7, test_trivial=false; end

% display something
if verbose, disp('  Doing prep work for emergency control');end
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
branch_st(isnan(branch_st)) = 1; % if the comm system has made this unreadable, assume branch is closed
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
nodes = (1:n)';
links = [F,T];
[grid_no,n_sub] = findSubGraphs(nodes,links);
% choose a bus that is already closest to zero
theta = ps.bus(:,C.bu.Vang) * pi / 180;
for grid_i = 1:n_sub
    bus_subset = find(grid_no == grid_i);
    theta_sub = theta(bus_subset);
    % find the bus in this island with the smallest absolute angle
    [~,ref_tmp] = min(abs(theta_sub));
    ref_ix = bus_subset(ref_tmp);
    % set this bus to be the local reference
    x_min(ix.x.d_theta(ref_ix)) = 0;
    x_max(ix.x.d_theta(ref_ix)) = 0;
end
cost = zeros(ix.nx,1);
cost(ix.x.dPd) = -1;
cost(ix.x.f_over) = f_over_cost;

%% set up equality constraints
% DC power flow constraint
A_pf = sparse(n,ix.nx); 
% generator injection constraint:
A_pf = A_pf + sparse(G,ix.x.dPg,1,n,ix.nx);
% load injection constraint:
A_pf = A_pf + sparse(D,ix.x.dPd,-1,n,ix.nx);
% DC power flow matrix constraint
% -B = sparse(F,T,-inv_X,n,n) + ...
%     sparse(T,F,-inv_X,n,n) + ...
%     sparse(T,T,+inv_X,n,n) + ...
%     sparse(F,F,+inv_X,n,n);
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
if verbose
    disp('  Solving the emergency control problem');
end
[x_star,~,exitflag,~] = cplexlp(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);

%% check and process the solution
if exitflag==1
    delta_Pd_pu = x_star(ix.x.dPd);
    delta_Pg_pu = x_star(ix.x.dPg);
    delta_Pg = delta_Pg_pu * ps.baseMVA;
    delta_Pd = delta_Pd_pu * ps.baseMVA;
    if verbose
        disp('  Solved the emergency control problem');
    end
else
    keyboard
    if verbose
        disp('Optimization failed');
        keyboard
    end
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
if any( Pd0>0 & ( Pd < 0 | Pd > Pd0 ) )
    error('Control results are outside of load limits');
    %p = find( Pd0>0 & ( Pd < 0 | Pd > Pd0 ) )
end
%flow_pu = measured_flow_pu + A_flow * x_star;
%keyboard
