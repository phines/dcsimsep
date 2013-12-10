clear all;

% ps = case30_ps;
ps = case6_ps;

C = psconstants;


% ps.branch([4 7 8 11],C.br.status) = 0;
ps.branch(:,C.br.status) = randi(2,1,size(ps.branch,1),1)'-1;

% ps.branch([1 3 5 10],C.br.status) = 0;

ps = dcpf(ps)

printps(ps)


% ps.shunt(23,1) = 9;
% ps.shunt(24,1) = 28; 

% ps.branch([4 10],C.br.Pf) = 80;
% ps.branch([15 24],C.br.Pf) = 70;
% ps.branch(:,C.br.rateA) = 10*randi(10,size(ps.branch,1),1);
ps.branch(:,C.br.rateA) = 80*rand(size(ps.branch,1),1);

% ps.branch(:,C.br.rates) = 10;


EPS = 1e-6;
br_st = ps.branch(:,C.br.status)~=0; % we could add conditions to this...

measured_flow = ps.branch(br_st(:),C.br.Pf);

% % check the inputs
% if nargin<2, error('need at least 2 inputs'); end
% if nargin<3, verbose = true; end
measured_flow_pu = measured_flow / ps.baseMVA;
% if verbose
%     disp('Doing prep work');
% end

% number of buses in the system
n = size(ps.bus,1);

% collect load data
D = ps.bus(ps.shunt(:,1)); % load bus locations
nd = length(D);
Pd0 = ps.shunt(:,C.sh.P) ./ ps.baseMVA; % in per unit
% collect generator data
G = ps.bus(ps.gen(:,1));  % gen bus locations
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
F = full(ps.bus(ps.branch(br_st,1)));
T = full(ps.bus(ps.branch(br_st,2)));
X = ps.branch(br_st,C.br.X);
inv_X = 1./X;
flow_max = ps.branch(br_st,C.br.rateA) / ps.baseMVA; % in per unit
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
% if verbose
%     disp('Solving the problem');
% end
% x_star = linprog(cost,A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);
x_star = cplexlp(cost',A_ineq,b_ineq,A_pf,b_pf,x_min,x_max);


delta_Pd_pu = x_star(ix.x.dPd);
delta_Pg_pu = x_star(ix.x.dPg);
delta_Pg = delta_Pg_pu * ps.baseMVA;
delta_Pd = delta_Pd_pu * ps.baseMVA;

ps.gen(:,C.gen.Pg) = ps.gen(:,C.gen.Pg) + delta_Pg;

ps.shunt(:,C.shunt.P) = ps.shunt(:,C.shunt.P) + delta_Pd;


ps = dcpf(ps)
printps(ps)
