
clear all;

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

br_outages = [22, 24, 27, 30];

opt = psoptions;
opt.verbose=true;

ps = ps_50;
Pdsum = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.status))
Pgsum = sum(ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status))


[is_blackout,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);


fprintf('\n\n\n')

ps = ps_119;
[is_blackout,relay_outages,MW_lost,p_out,busessep,flows] = dcsimsep(ps,br_outages,[],opt);

