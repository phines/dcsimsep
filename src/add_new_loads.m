function shunt = add_new_loads(ps)
% adds new loads so that there is one load at each bus

error('Don''t use this file');

n = size(ps.bus,1);
C = psconstants;

% find the buses with no loads
has_load = unique(ps.bus_i(ps.shunt(:,1)));
does_not_have_load = setdiff(1:n,has_load);

% create new loads at these locations
n_new = length(does_not_have_load);
new_shunts = zeros(n_new,C.sh.cols);
bus_nos = ps.bus(does_not_have_load,1);
new_shunts(:,C.sh.bus) = bus_nos;
new_shunts(:,C.sh.status) = 1;
new_shunts(:,C.sh.frac_S) = 1;
shunt = [ps.shunt;new_shunts];