function shunt_index = assign_loads_to_buses(ps)
% usage: shunt_index = assign_loads_to_buses(ps)

n = size(ps.bus,1);
nd = size(ps.shunt,1);

shunt_index = nan(n,1);

% assign all of the obvious ones
for d = 1:nd
    bi = ps.bus_i(ps.shunt(:,1));
    shunt_index(bi) = d;
end

% asign all of the others
% some other day...
