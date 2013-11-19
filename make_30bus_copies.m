% make copies of 30 bus case case

ps0 = case30_mod_ps;
ps0 = dcpf(ps0);
ps1 = ps0;
ps2 = ps0;
ps3 = ps0;
ps4 = ps0;
% fix bus numbers
ps1.bus(:,1) = ps1.bus(:,1)+100;
ps2.bus(:,1) = ps2.bus(:,1)+200;
ps3.bus(:,1) = ps3.bus(:,1)+300;
ps4.bus(:,1) = ps4.bus(:,1)+400;

% fix the locations
ps2.bus(:,C.bu.locX) = ps1.bus(:,C.bu.locX) + 1;
ps3.bus(:,C.bu.locY) = ps1.bus(:,C.bu.locY) + 1;
ps4.bus(:,C.bu.locX) = ps1.bus(:,C.bu.locX) + 1;
ps4.bus(:,C.bu.locY) = ps1.bus(:,C.bu.locY) + 1;
% merge the cases
ps = ps0;
ps.bus = [ps1.bus;bs2.bus;ps3.bus;ps4.bus];
ps.bus = [ps1.bus;bs2.bus;ps3.bus;ps4.bus];

% add new branches to link cases
% from 1 - 2
% from 1 - 3
% from 3 - 4
% from 2 - 4

% 