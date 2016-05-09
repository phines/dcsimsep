function print_load_shedding_results(ps,delta_Pg,delta_Pd)
C = psconstants;
EPS = 1e-4;
% generation shedding
fprintf('   Generation shedding (%.4g MW total):\n',-sum(delta_Pg));
ramp_set = find(delta_Pg < -EPS);
all_ramp_bus = unique(ps.gen(ramp_set,C.ge.bus));
for j = 1:length(all_ramp_bus)
    idx_ramp_set = ps.gen(ramp_set,C.ge.bus) == all_ramp_bus(j);
    idx_Pg = ramp_set(idx_ramp_set);
    this_ramp = -sum(delta_Pg(idx_Pg));
    fprintf('    %.4g MW on bus %d.\n',this_ramp,all_ramp_bus(j));
end
% load shedding
fprintf('   Load shedding (%.4g MW total):\n',-sum(delta_Pd));
shed_set = find(delta_Pd < -EPS);
all_shed_bus = unique(ps.shunt(shed_set,C.sh.bus));
for j = 1:length(all_shed_bus)
    idx_shed_set = ps.shunt(shed_set,C.ge.bus) == all_shed_bus(j);
    idx_Pd = shed_set(idx_shed_set);
    this_shed = -sum(delta_Pd(idx_Pd));
    fprintf('    %.4g MW on bus %d.\n',this_shed,all_shed_bus(j));
end
 