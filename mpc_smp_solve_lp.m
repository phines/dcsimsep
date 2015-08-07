function [dPd, dPg] = mpc_smp_solve_lp(ps_agent,OptVarBus_I,ramp_limits,opt) 
% solves an MPC model for the stress mitigation problem (SMP) at each 
% sample using solvelp. Generation and load can only vary on the local 
% neighborhood. In addition, only branches connected to local buses are 
% constrained by fmax.

C = psconstants;
% get some information
Np = opt.sim.Np;
% collect some data
nbus = size(ps_agent.bus,1);
br_st = (ps_agent.branch(:,C.br.status) == 1);
F = full(ps_agent.bus_i(ps_agent.branch(br_st,1)));
T = full(ps_agent.bus_i(ps_agent.branch(br_st,2)));
X = ps_agent.branch(br_st,C.br.X);
G = ps_agent.bus_i(ps_agent.gen(:,1));    % generator locations
D = ps_agent.bus_i(ps_agent.shunt(:,1));  % demand/shunt locations
% find all the generators, loads and branches in the local neighborhood
if isempty(OptVarBus_I)
    OptVarBus_I = 1:nbus;
end
ig = ismember(G, OptVarBus_I);
ng = sum(ig);
ish = ismember(D, OptVarBus_I);
nd = sum(ish);
if ng == 0 || nd == 0
    dPd = 0; 
    dPg = 0;
    return
end
Pg0 = ps_agent.gen(ig,C.ge.P) .* ps_agent.gen(ig,C.ge.status);
Pd0 = ps_agent.shunt(ish,C.sh.P) .* ps_agent.shunt(ish,C.sh.factor);
[br_id,~] = find(ismember(F,OptVarBus_I) | ismember(T,OptVarBus_I));
nbr = length(br_id);
f0_all = ps_agent.branch(br_st,C.br.Pf)/ps_agent.baseMVA;
f0 = f0_all(br_id);
fmax_all = ps_agent.branch(br_st,C.br.rateB)/ps_agent.baseMVA;
fmax = fmax_all(br_id);
% assign one ref bus for each island: choose a bus that is already closest to zero
nodes = (1:nbus)'; links = [F,T];
grid_no = FindSubGraphs(nodes,links);
n_sub = max(grid_no);
ref_bus_i = zeros(n_sub,1);
theta = ps_agent.bus(:,C.bu.Vang) * pi / 180;
for grid_i = 1:n_sub
    bus_subset = find(grid_no == grid_i);
    theta_sub = theta(bus_subset);
    % find the bus in this island with the smallest absolute angle
    [~,sub_idx] = min(abs(theta_sub));
    % Add this bus to the reference list
    ref_bus_i(grid_i) = bus_subset(sub_idx);
end
% build an index into the decesion vector
ix.dPg    = reshape(1:ng*Np, ng, Np);
ix.dPd    = reshape(1:nd*Np, nd, Np)     + ng * Np;
ix.dtheta = reshape(1:nbus*Np, nbus, Np) + (ng + nd) * Np;
ix.s      = reshape(1:nbr*Np, nbr, Np)   + (ng + nd + nbus) * Np;
ix.f      = reshape(1:nbr*Np, nbr, Np)   + (ng + nd + nbus + nbr) * Np;
nx = (ng + nd + nbus + nbr + nbr) * Np;

%% set up x_min, x_max
x_min = nan(nx,1);
x_max = nan(nx,1);
% dPg
x_max(ix.dPg) = 0;
x_min(ix.dPg) = [-ps_agent.gen(ig,C.ge.ramp_rate_down); repmat(-ramp_limits(ig),Np-1,1)]; % ramping for the first time step is limited based on the agent information 
% dPd
x_max(ix.dPd) = 0;
x_min(ix.dPd) = -Inf;
% dtheta
x_min(ix.dtheta) = -Inf; x_max(ix.dtheta) = Inf;
x_min(ix.dtheta(ref_bus_i,:)) = 0;
x_max(ix.dtheta(ref_bus_i,:)) = 0;
% s
x_min(ix.s) = 0;
x_max(ix.s) = Inf;
% f
x_min(ix.f) = -Inf;
x_max(ix.f) = Inf;

%% set up the objective
c_obj = zeros(nx,1);
c_obj(ix.dPd) = -opt.sim.cost.load / ps_agent.baseMVA; % changed to pu in obj func, just to make it compatible with emergency_control_dec
c_obj(ix.s) = opt.sim.cost.overload;

%% Constraints
% dc power flow constraints
inv_X = (1./X);
B = sparse(F,T,-inv_X,nbus,nbus) + ...
    sparse(T,F,-inv_X,nbus,nbus) + ...
    sparse(T,T,+inv_X,nbus,nbus) + ...
    sparse(F,F,+inv_X,nbus,nbus);
extra_col = ix.dtheta(1) - 1;
B_bd = make_block_diag(B,Np,extra_col,nx);
Ag = sparse(G(ig),1:ng,1,nbus,ng)/ps_agent.baseMVA; % To make dPg and dPd in pu
extra_col = ix.dPg(1) - 1;
Ag_bd = make_block_diag(Ag,Np,extra_col,nx);
Ad = sparse(D(ish),1:nd,1,nbus,nd)/ps_agent.baseMVA;
extra_col = ix.dPd(1) - 1;
Ad_bd = make_block_diag(Ad,Np,extra_col,nx);

A_dc = B_bd - Ag_bd + Ad_bd;
b_dc = sparse(nbus*Np,1);
% update b_dc for the first time step based on the agent information ->
% look at emergency_control_dec for more information
dPg = ps_agent.delta_Pg0_all / ps_agent.baseMVA;
dPd = ps_agent.delta_Pd0_all / ps_agent.baseMVA;
b_dc_1 = sparse(G,1,-dPg,nbus,1) + sparse(D,1,dPd,nbus,1);
mismatch = sum(dPg) - sum(dPd); 
% make b_pf balanced when it is not already so (add mismatch to load/gen on the agent itself)
this_bus_i = ps_agent.bus_i(ps_agent.bus_id);
b_dc_1(this_bus_i) = b_dc_1(this_bus_i) + mismatch;
% update b_dc
b_dc(1:nbus) = -b_dc_1;

% flow computation constraints
F_sub = F(br_id);
T_sub = T(br_id);
X_sub = X(br_id);
A_f = sparse(1:nbr*Np, ix.f, 1, nbr*Np, nx) + ...
      sparse(nbr+1:nbr*Np, ix.f(:,1:end-1), -1, nbr*Np, nx) + ...
      sparse(1:nbr*Np, ix.dtheta(F_sub,:), repmat(-1./X_sub,1,Np), nbr*Np, nx) + ...
      sparse(1:nbr*Np, ix.dtheta(T_sub,:), repmat(1./X_sub,1,Np), nbr*Np, nx);

b_f = [f0; zeros(nbr*(Np-1),1)];

% flow - s <= fmax
A_flim1 = sparse(1:nbr*Np, ix.f, 1, nbr*Np, nx) + ...
          sparse(1:nbr*Np, ix.s, -1, nbr*Np, nx);
b_flim1 = repmat(fmax,Np,1);

% -flow - s <= fmax
A_flim2 = sparse(1:nbr*Np, ix.f, -1, nbr*Np, nx) + ...
          sparse(1:nbr*Np, ix.s, -1, nbr*Np, nx);
b_flim2 = repmat(fmax,Np,1);

% -\sigma dPg <= Pg0
A_dPg = sparse(repmat((1:ng)',Np,1), ix.dPg, -1, ng, nx);
b_dPg = Pg0;

% -\sigma dPd <= Pd0
A_dPd = sparse(repmat((1:nd)',Np,1), ix.dPd, -1, nd, nx);
b_dPd = Pd0;

% merge the constraints
lp_Aeq = [A_dc; A_f];
lp_beq = [b_dc; b_f];
lp_Aineq = [A_flim1; A_flim2; A_dPg; A_dPd];
lp_bineq = [b_flim1; b_flim2; b_dPg; b_dPd];

%% solve the problem
options = cplexoptimset;
[x,fval,exitflag,output] = ...
    cplexlp(c_obj,lp_Aineq,lp_bineq,lp_Aeq,lp_beq,x_min,x_max,[],options);
if exitflag > 0
    dPg = zeros(size(ps_agent.gen(:,1))); 
    dPd = zeros(size(ps_agent.shunt(:,1)));
    dPg(ig) = x(ix.dPg(1:ng));
    dPd(ish) = x(ix.dPd(1:nd));
    if opt.verbose
        fprintf('  Solved the emergency control problem on bus %d using MPC. fval = %g\n', ...
            ps_agent.bus_id, fval);
    end
end




    

