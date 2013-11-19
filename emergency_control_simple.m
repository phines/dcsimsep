function [delta_Pd,delta_Pg] = emergency_control_simple(ps,measured_flow,verbose)
% An emergency control function to adjust load, based on current measured
%  flows
% usage: [delta_df,delta_Pg] = emergency_opf(ps,measured_flow,verbose)
% inputs:
%  ps - power systems data
%  measured_flow - the amount of flow (in MW/MVA) on each transmission line
%  verbose - a binary indicating whether to print stuff about what is going
%   on.
%  
% outputs:
%  delta_Pd - changes to amount of load at each "shunt"
%  delta_Pg - changes to the generator output at each generator

%% Prep work
C = psconstants;
EPS = 1e-6;

% check the inputs
if nargin<2, error('need at least 2 inputs'); end
if nargin<3, verbose = true; end
measured_flow_pu = measured_flow / ps.baseMVA;
if verbose
    disp('Doing prep work');
end

% number of buses in the system
n = size(ps.bus,1);

% collect load data
D = ps.bus_i(ps.shunt(:,1)); % load bus locations
nd = length(D);
Pd0 = ps.shunt(:,C.sh.P) ./ ps.baseMVA; % in per unit
% collect generator data
G = ps.bus_i(ps.gen(:,1));  % gen bus locations
ng = length(G);
ge_status = ps.gen(:,C.ge.status)==1;
Pg0 = ps.gen(:,C.ge.P).*ge_status./ ps.baseMVA; % in per unit
Pg = Pg0;
Pmax = Pg0; % < this could change
Pmin = zeros(ng,1);
if any(Pg<Pmin-EPS) || any(Pg>Pmax+EPS)
    error('Generation outside of Pmin/Pmax');
end
% collect transmission line data
br_st = ps.branch(:,C.br.status)~=0; % we could add conditions to this...
F = full(ps.bus_i(ps.branch(br_st,1)));
T = full(ps.bus_i(ps.branch(br_st,2)));
X = ps.branch(br_st,C.br.X);
inv_X = 1./X;
flow_max = ps.branch(br_st,C.br.rateB) / ps.baseMVA; % in per unit
m = length(F); % number of transmission lines in the system

%% set up an index so that we can find things
ix.x.d_theta = (1:n);
ix.x.dPg = (1:ng) + n;
ix.x.dPd = (1:nd) + n + ng;
ix.nx = n + ng + nd;
% we might also need an index for the constraints...???

%% set up limits and costs for the variables
x_min = zeros(ix.nx,1)-Inf;
x_max = zeros(ix.nx,1)+Inf;
% for dPg
x_min(ix.x.dPg) = -Pg0;
x_max(ix.x.dPg) = 0;
% for dPd
x_min(ix.x.dPd) = -Pd0;
x_max(ix.x.dPd) = 0;
% for theta
x_min(ix.x.d_theta(1)) = 0;
x_max(ix.x.d_theta(1)) = 0;
cost = zeros(ix.nx,1);
cost(ix.x.dPd) = -1;

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
A_flow = sparse((1:m),F,inv_X,m,ix.nx) + ...
         sparse((1:m),T,-inv_X,m,ix.nx);
b_flow_left  = -flow_max - measured_flow_pu;
b_flow_right = +flow_max - measured_flow_pu;
A_ineq = [-A_flow;A_flow];
b_ineq = [-b_flow_left;b_flow_right];

%% run the optimization
if verbose
    disp('Solving the problem');
end
x_star = linprog(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);

%% check and process the solution
delta_Pd_pu = x_star(ix.x.dPd);
delta_Pg_pu = x_star(ix.x.dPg);
delta_Pg = delta_Pg_pu * ps.baseMVA;
delta_Pd = delta_Pd_pu * ps.baseMVA;

%flow_pu = measured_flow_pu + A_flow * x_star;
%keyboard
