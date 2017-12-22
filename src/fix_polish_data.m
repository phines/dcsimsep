% fix the polish data
close all;
C = psconstants;
ps = updateps(case2383wp_ps);
ps.branch(:,C.br.status) = 1;
ps.shunt(:,C.sh.P) = ps.shunt(:,C.sh.P) * 1.1;
n = size(ps.bus,1);
ramp_rate = ps.gen(:,C.ge.Pmax);

% run the initial dcpf
C = psconstants;
ps = rebalance(ps,ones(n,1),ramp_rate);
ps = dcpf(ps);
ps_base = ps; % save the base case

rateA = ps.branch(:,C.br.rateA);
rateB = ps.branch(:,C.br.rateB);
rateC = ps.branch(:,C.br.rateC);
flows = max(abs(ps.branch(:,[C.br.Pf C.br.Pt])),[],2);
figure(1); clf;
load_ratio = flows./rateA;
hist(load_ratio,50);

rateA_ = max(ceil(flows*1.05),rateA);

figure(2); clf;
hist(flows./rateA_,50);
title('base case over rate A');
return

% run all single branch contingencies
contingency_flows = flows;
m = size(ps.branch,1);
ng = size(ps.gen,1);
for i=1:m
    fprintf('%d of %d\n',i,m);
    % return to the base case
    ps = ps_base;
    % take out the transmission line
    ps.branch(i,C.br.status) = 0;
    % find the subgrids
    br_st = ps.branch(:,C.br.status)~=0;
    [sub_grids,n_sub] = find_subgraphs(ps.bus(:,1),ps.branch(br_st,1:2));
    % rebalance if needed
    if n_sub>1
        ps = rebalance(ps,sub_grids,ramp_rate);
    end
    % run the dcpf
    ps = dcpf(ps,sub_grids,false);
    flow = max(abs(ps.branch(:,[C.br.Pf C.br.Pt])),[],2);
    contingency_flows = max(flow,contingency_flows);
end
% edit rateB
rateB_ = max(ceil(contingency_flows*1.01),rateB);
figure(3); clf;
hist(contingency_flows./rateB_,50);
title('contingency over rate B');

rateC_ = max(rateB_*(1.5/1.2),rateC);
figure(3); clf;
hist(contingency_flows./rateC_,50);
title('contingency over rate C');

%% now fix the data

% revert to the base case
ps = ps_base;
% record the revised limits
ps.branch(:,C.br.rateA) = rateA_;
ps.branch(:,C.br.rateB) = rateB_;
ps.branch(:,C.br.rateC) = rateC_;

% write out the case data
ps.branch(:,C.br.status) = 1;
writeps(ps,'case2383_mod_ps');

save temp;

%% look at the case data
load temp ;
figure(1);
hist(rateA_./rateA,100);
figure(2);
hist(rateB_./rateB,100);

return

%% validate the case
clear all
C = psconstants;
ps = case2383_mod_ps;
ps = dcpf(ps);
ramp_rate = ps.gen(:,C.ge.Pmax);
% run all single branch contingencies
m = size(ps.branch,1);
contingency_flows = zeros(m,1);
ng = size(ps.gen,1);
rateB = ps.branch(:,C.br.rateB);

for i=1:m
    % take out the transmission line
    ps.branch(:,C.br.status) = 1;
    ps.branch(i,C.br.status) = 0;
    % find the subgrids
    br_st = ps.branch(:,C.br.status)~=0;
    [sub_grids,n_sub] = find_subgraphs(ps.bus(:,1),ps.branch(br_st,1:2));
    % rebalance if needed
    if n_sub>1
        ps = rebalance(ps,sub_grids,ramp_rate);
    end
    % run the dclf
    ps = dcpf(ps,sub_grids);
    flow = max(abs(ps.branch(:,[C.br.Pf C.br.Pt])),[],2);
    if any(flow>rateB)
        fprintf('Branch %d resulted in a violation\n',i);
        keyboard
    end
    fprintf('Checked %d of %d\n',i,m);
end

