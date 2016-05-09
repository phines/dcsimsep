function is_powered = find_buses_with_power(ps,opt)
% figure out which buses in the system have power to them

C = psconstants;

% initialize the ouptut
is_powered = ps.bus(:,C.bu.status);
n = length(is_powered);

% look at each bus
for i = 1:n
    % if we we are doing 2-way coupling, and this bus has power
    if opt.comm.two_way && is_powered(i)
        % figure out which bus this bus is powered from
        sh_index = ps.bus(i,C.bu.power_from_sh);
        p_failure = 1-ps.shunt(sh_index,C.sh.factor);
        is_failed = rand<p_failure;
        if is_failed
            is_powered(i) = false;
        end
    end
end

return

%{ 
OLD STUFF
load_shedding_limit = 0.5; % < This is a big assumption!
n = size(ps.bus,1);
is_powered = true(n,1);
G = ps.bus_i(ps.gen(:,1));
ge_status = ps.gen(:,C.ge.status)==1;
D = ps.bus_i(ps.shunt(:,1));
d_factor = ps.shunt(:,C.sh.factor);
Pg0 = ps.gen(:,C.ge.P).*ge_status;
Pg  = Pg0;
EPS = 1e-6;

% go through each subgrid
grid_list = unique(sub_grids)';
for g = grid_list
    bus_set = find(g==sub_grids);
    % find the amount of generation:
    Gsub = ismember(G,bus_set) & ge_status;
    Pg_sub = sum(Pg(Gsub));
    % find the amount of load shedding:
    Dsub = ismember(D,bus_set);
    mean_load_shedding = mean(1-d_factor(Dsub));
    % if there is no active generation in this subgrid, or there is a lot of load shedding, then no power
    if Pg_sub<EPS || mean_load_shedding > load_shedding_limit
        is_powered(bus_set)=false;
    end
end
%}