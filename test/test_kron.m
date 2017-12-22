
ps = case30_ps;
C = psconstants;

ps = rebalance(ps);
ps = dcpf(ps);
printps(ps);

%% Code extracted from dcpf
% extract some data
n = size(ps.bus,1);
m = size(ps.branch,1);
br_st = ps.branch(:,C.br.status)~=0;
F = full(ps.bus_i(ps.branch(br_st,1)));
T = full(ps.bus_i(ps.branch(br_st,2)));
X = ps.branch(br_st,C.br.X);
G = ps.bus_i(ps.gen(:,1));    % generator locations
D = ps.bus_i(ps.shunt(:,1));  % demand/shunt locations
sf = ps.shunt(:,C.sh.status); % shunt factor, used for load shedding
Pg_pu = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status) / ps.baseMVA;
Pg_max_pu = ps.gen(:,C.ge.Pmax).*ps.gen(:,C.ge.status) / ps.baseMVA;
Pd_pu = ps.shunt(:,C.sh.P) / ps.baseMVA;

% calculate the B matrix and initial theta
theta = zeros(n,1);
Vmag = ones(n,1);
inv_X = (1./X);
B = sparse(F,T,-inv_X,n,n) + ...
    sparse(T,F,-inv_X,n,n) + ...
    sparse(T,T,+inv_X,n,n) + ...
    sparse(F,F,+inv_X,n,n);

% find the net generation
Pg_max_full = full(sparse(G,1,Pg_max_pu,n,1));
Pg_full = full(sparse(G,1,Pg_pu,n,1));
Pg_org  = Pg_full;
Pd_full = full(sparse(D,1,Pd_pu.*sf,n,1));
net_gen = Pg_full - Pd_full;

% check the dcpf from here
ref = find(ps.bus(:,C.bu.type)==C.REF);
non_ref = setdiff(1:n,ref);
theta(non_ref) = B(non_ref,non_ref)\net_gen(non_ref);

%% Do the Kron reduction
K = [2:5:30]
U = setdiff(1:30,K)

Bkron = B(K,K) - ( B(K,U) * inv(B(U,U)) ) * B(U,K)

%% test the Kron reduction
% increase gen at the first keep bus and increase load at the last one
d_net_gen = zeros(n,1);
d_net_gen(K(1)) = 0.01;
d_net_gen(K(end)) = -0.01;
net_gen = net_gen + d_net_gen;

theta_new = zeros(n,1);
theta_new(non_ref) = B(non_ref,non_ref)\net_gen(non_ref);

d_theta_kron = zeros(n,1);
d_theta_kron(K) = Bkron \ d_net_gen(K)
d_theta = theta_new - theta;

[theta theta_new d_theta d_theta_kron]