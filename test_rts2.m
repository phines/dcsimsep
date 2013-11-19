
clear all;
clc;

% ps = mp2ps(case73_rts_96_modified);

% ps = dcopf(ps);
% 
% m = size(ps.branch,1);
% 
% for i=1:m
%     ps.branch(:,C.br.status) = 1;
%     ps.branch(i,C.br.status) = 0;
%     ps = dcpf(ps);
%     printps(ps);
%     pause
% end

% return
load ps_RTS_all
C = psconstants;
ps0 = mp2ps(case73_rts_96_modified);
% P0 = ps0.shunt(:,C.sh.P);
% P1 = ps_50.shunt(:,C.sh.P);

br_outages = [22, 24, 27, 30, 28];

opt = psoptions;
opt.verbose=true;


ps = ps_50;
% trip branches
ps.branch(br_outages,C.br.status) = 0;
ps = dcpf(ps);
flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
Pg = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status);
Pg_max = ps.gen(:,C.ge.Pmax).*ps.gen(:,C.ge.status);
figure(2);
hist(abs(flow)./flow_max,100);
Pd0 = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.status))
title('50');
return
fprintf('\n\n\n')


ps = ps_119;
ps.branch(br_outages,C.br.status) = 0;
ps = dcpf(ps);
flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
Pg = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status);
Pg_max = ps.gen(:,C.ge.Pmax).*ps.gen(:,C.ge.status);
figure(1);
hist(abs(flow)./flow_max,100);
title('119');
Pd1 = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.status))

