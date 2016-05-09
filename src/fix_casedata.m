%casename = 'case300_001_ps');
casename = 'case2383_mod_ps';

% extract data
C  = psconstants;
ps = updateps(feval(casename));
ps.branch(:,C.br.status) = 1;
ps = dcpf(ps);
n = size(ps.bus,1);
m = size(ps.branch,1);
worst_flow = zeros(m,1);

% find the worst case single contingency flows
flow_limit
for br = 1:m
    % switch off the branch
    ps.branch(i,br.status) = 0;
    % run the power flow
    ps = dcpf(ps);
    % check the branch flows
    flow = abs(ps.branch(:,C.br.Pf));
    % update the worst case flow
    worst_flow = max(flow,worst_flow);
    % display something
    fprintf('Completed branch %d of %d\n',
end
