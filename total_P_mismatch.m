function mis = total_P_mismatch(ps)
% a simple function that computes the total system mismatch
C = psconstants;
Pd = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
Pg = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status);

mis = sum(Pg) - sum(Pd);
