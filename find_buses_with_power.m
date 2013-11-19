function is_powered = find_buses_with_power(ps,sub_grids)
% figure out which buses in the system have power to them

% initialize the ouptut
C = psconstants;
n = size(ps.bus,1);
is_powered = true(n,1);
G = ps.bus_i(ps.gen(:,1));
ge_status = ps.gen(:,C.ge.status)==1;
Pg0 = ps.gen(:,C.ge.P).*ge_status;
Pg  = Pg0;
EPS = 1e-6;

% go through each subgrid
grid_list = unique(sub_grids)';
for g = grid_list
    bus_set = find(g==sub_grids);
    Gsub = ismember(G,bus_set) & ge_status;
    %Dsub = ismember(D,bus_set) & d_factor>0;
    Pg_sub = sum(Pg(Gsub));
    % if there is no active generation in this subgrid, then no power
    if Pg_sub<EPS
        is_powered(bus_set)=false;
    end
end
