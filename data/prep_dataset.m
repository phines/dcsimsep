% Use this script to prepare a dataset for use with dcsimsep

% load the initial dataset
C = psconstants;
ps = updateps(case2383_mod_ps);
ps.bus(:,C.bu.comm_status) = 1;
ps.bus(:,C.bu.status) = 1;
ps.bus(:,C.bu.power_from_sh) = assign_loads_to_buses(ps);
