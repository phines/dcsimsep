function [ps,sub_grids,n_sub] = dcpf(ps,sub_grids,load_shedding)
% usage: [ps,sub_grids,n_sub] = dcpf(ps,sub_grids,load_shedding)
% a very simple dc power flow calculation
% if sub_grids are not specified, they are calculated from the graph.

% input check
if nargin<1, error('ps structure must be specified'); end;
if nargin<2, sub_grids = []; end
if nargin<3, load_shedding=false; end
verbose = false;

% some constants
C = psconstants;
ps = updateps(ps);
EPS = 1e-6;

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

% find the ref bus
ref = find(ps.bus(:,C.bu.type)==C.REF);
if length(ref)~=1 || Pg_full(ref)<=0
    [~,ref] = max(Pg_full);
end

% find sub-grids
if isempty(sub_grids)
    [sub_grids,n_sub] = findSubGraphs(ps.bus(:,1),ps.branch(br_st,1:2));
else
    n_sub = length(unique(sub_grids));
end
ref_record = zeros(1,n_sub); % to save the reference buses

for g = 1:n_sub
    % figure out which buses are in this subgrid
    if n_sub>1
        subset = (sub_grids==g);
    else
        subset = true(n,1);
    end
    if ~any(subset), error('no buses in this subset???'); end
    % find a reference bus
    if ~subset(ref)
        % find the largest generator in the subset
        [~,ref] = max(Pg_full.*subset);
    end
    if length(ref)>1
        error('multiple ref buses');
    end
    % check for a blackout in this subgrid
    if sum(Pd_full(subset))<=0 || sum(Pg_max_full(subset))<=0 %|| Pg_full(ref)<=0
        if load_shedding
            theta(subset) = 0;
            Pg_full(subset) = 0;
            Pd_full(subset) = 0;
            Vmag(subset) = 0;
            net_gen(subset) = 0;
            % shut off shunts
            bus_list = find(subset);
            sh_subset = ismember(D,bus_list);
            sf(sh_subset) = 0;
            % shut off generators
            ge_subset = ismember(G,bus_list);
            ps.gen(ge_subset,C.ge.status) = 0;
            if verbose
                fprintf('dcpf found a blackout in subgrid %d of %d\n',g, n_sub);
            end
        else
            error('Power flow failed to converge');
        end
        continue
    end
    % remove the slack bus from the subset in order to find angles
    subset(ref) = false;
    % find the angles
    theta(subset) = B(subset,subset)\net_gen(subset);
    % record the reference bus
    ref_record(g) = ref;
end
% find the new power flows
Pft_pu = zeros(m,1);
Pft_pu(br_st) = inv_X .* (theta(F) - theta(T));
% record the flows to the ps structure
ps.branch(:,C.br.Imag_f) = abs(Pft_pu);
ps.branch(:,C.br.Imag_t) = abs(Pft_pu);
ps.branch(:,C.br.Pf) = +Pft_pu * ps.baseMVA;
ps.branch(:,C.br.Pt) = -Pft_pu * ps.baseMVA;
ps.branch(:,C.br.Qf) = 0; % dcpf assumption
ps.branch(:,C.br.Qt) = 0; % dcpf assumption
% record the voltages to the ps structure
ps.bus(:,C.bu.Vang) = theta*180/pi;
ps.bus(:,C.bu.Vmag) = Vmag;
% fix the generation imbalance, if any
Pg_full = B*theta + Pd_full;
Pg_full = round(Pg_full*1e6)/1e6;
delta_Pg = Pg_full - Pg_org;
subset = find(abs(delta_Pg)>EPS)';
for b = subset
    % find the generators at this bus
    bus_no = ps.bus(b,1);
    genset = (ps.gen(:,1)==bus_no);
    if ~any(genset), error('no generators?'); end;
    % find the new Pg for these generators
    Pg_old = Pg_pu(genset);
    Pg_new = Pg_old * (Pg_full(b)/sum(Pg_old));
    % record into the generator data
    ps.gen(genset,C.ge.Pg) = Pg_new * ps.baseMVA;
end
ps.gen(:,C.ge.Qg) = 0; % dcpf assumption...
% record the load factor
ps.shunt(:,C.sh.factor) = sf;


% record the B matrix
ps.B = B;

% check for imbalance
imbalance = sum(ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status)) - sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
if abs(imbalance)>1e-3
    error('Imbalance found when there shouldn''t be one');
end

%{
% debug
if 0
    Imag = max(ps.branch(:,C.br.Imag_f),ps.branch(:,C.br.Imag_t));
    Imax = ps.branch(:,C.br.rateB)/ps.baseMVA;
    figure(1);
    hist(Imag./Imax,30);
    keyboard
end
%}
