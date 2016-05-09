function mis = total_P_mismatch(ps,sub_grids)
% a simple function that computes the total system mismatch in MW

% constants
C = psconstants;
EPS = 1e-3;
% init the output
mis = 0;

% find the load and generation
Pd = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
Pg = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status);
mis = sum(Pg) - sum(Pd);
% if abs(mis)>EPS
%     keyboard
% end

if nargin>1
    % map the mismatch to buses
    n = size(ps.bus,1);
    G = ps.bus_i(ps.gen(:,1));    % generator locations
    D = ps.bus_i(ps.shunt(:,1));  % demand/shunt locations
    Pg_full = full(sparse(G,1,Pg,n,1));
    Pd_full = full(sparse(D,1,Pd,n,1));
    % figure out how many subgrids
    n_sub = length(unique(sub_grids));
    for g = 1:n_sub
        subset = (sub_grids==g);
        Pd_total = sum(Pd_full(subset));
        Pg_total = sum(Pg_full(subset));
        mis = max(abs(mis),abs(Pd_total - Pg_total));
        %if abs(mis)>EPS
        %    error('debug');
        %end
    end
end
