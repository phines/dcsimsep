function shunt_index = assign_loads_to_buses(ps)
% usage: shunt_index = assign_loads_to_buses(ps)

% prep work
C = psconstants;
n = size(ps.bus,1);
nd = size(ps.shunt,1);
% make the list of load buses
D = ps.bus_i(ps.shunt(:,1));

% initialize the output
shunt_index = nan(n,1);

% assign all of the obvious ones
for d = 1:nd
    % Find the bus that this load serves
    bi = ps.bus_i(ps.shunt(d,1));
    % Assign that bus to this load
    shunt_index(bi) = d;
end

% build the adjacency matrix
nodes = 1:n;
F = ps.bus_i(ps.branch(:,1));
T = ps.bus_i(ps.branch(:,2));
A = sparse([F;T],[T;F],1,n,n);

% Go through each bus that doesn't have a shunt, and give it one

bu_subset = find(isnan(shunt_index));
for i = bu_subset'
    potential_loads = [];
    % keep searching until we find a candidate load
    adj_buses = i;
    while isempty(potential_loads)
        % Find the buses located near this one
        [adj_buses,~] = find(A(:,adj_buses));
        if length(adj_buses)>1
            potential_loads = find(ismember(D,adj_buses));
            % if we found a shunt, then
            % choose a random shunt, apply it to this bus
            if ~isempty(potential_loads)
                ix = randi([1 length(potential_loads)]);
                shunt_index(i) = potential_loads(ix);
                break;
            end
        end
    end
end
