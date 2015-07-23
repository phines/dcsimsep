load case2383_mod_ps;
C = psconstants;
ps.gen(:,C.ge.Pg) = ps.gen(:,C.ge.Pg)*1.2;
ps.shunt(:,C.sh.P) = ps.shunt(:,C.sh.P)*1.2;
ps = dcpf(ps);
opt = psoptions;
opt.draw.width = 0.005;
opt.draw.bus_nos = false;
opt.draw.simple = true;
opt.draw.fontsize = 14;
figure(1); clf;
drawps(ps,opt);

