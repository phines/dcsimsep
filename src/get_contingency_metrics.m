function [metric,is_overload] = get_contingency_metrics(ps,branch_outages,method,LODF,weights)
% generate metrics that rank branch outage contingencies.
% usage: metric = get_contingency_metrics(ps,branch_outages,method,LODF)
% inputs:
%  ps - a power system structure
%  branch_outages - a (n_contingency) by (contingency size, n-1, n-2, etc) 
%   list of contingencies to evaluate.
%  method should be one of:
%   wollenberg - the method from Ejebe and Wollenberg, modified to use the LODF factors 
%    and n-k outages, using the method from Davis and Overbye
%   overbye - the method from Davis and Overbye
%   PCP - a "potential cascading power" metric devised by me.
%   max-ratio - the maximum of the flow/rateA ratios
%  ouputs:
%   metric
%   is_overload, indicates if each contingency results in an overload

% check the inputs
if nargin<3
    method = 'wollenberg';
end
if nargin<4
    % build the LODF matrix
    H = makePTDF(ps.baseMVA, ps.bus, ps.branch);
    LODF = makeLODF(ps.branch, H);
    LODF(isnan(LODF)) = 0;
    LODF(isinf(LODF)) = 0;
end
if nargin<5
    weights=1;
end

% extract data
C = psconstants;
n_outages = size(branch_outages,1);
k         = size(branch_outages,2);
m         = size(ps.branch,1);
% init the output
metric = zeros(n_outages,1);
is_overload = zeros(size(metric));
% extract the initial power flows
%  (assumes the power flow was performed externally)
flow = ps.branch(:,C.br.Pf)/ps.baseMVA;
flow_max = ps.branch(:,C.br.rateA)/ps.baseMVA;
if any(flow>flow_max+EPS)
    error('flows');
end
base_woll_metric = sum( weights.*(flow./flow_max).^2 );

% estimate the post contingency flows using Overbye's M
for i = 1:n_outages
    outage = branch_outages(i,:);
    if k==1
        flow_post = flow + LODF * sparse(outage,1,F,m,1);
    else
        M = -LODF(outage,outage);
        F = flow(outage);
        if cond(M)>1e6
            % these are cases where one of the outages is on a radial line
            flow_post = flow; % who knows what to do here...
        else
            flow_post = flow + LODF(:,outage) * (M\F);
        end
    end
    is_overload(i) = any(flow_post > flow_max);
    switch lower(method)
        case 'wollenberg'
            post_ratio = flow_post ./ flow_max;
            metric(i) = sum( weights.*(post_ratio.^2) ) - base_woll_metric;
        case 'pcp'
            is_in_excess = flow_post >= flow_max;
            metric(i) = sum( abs(flow_post) .* is_in_excess );
            % could also use a sigmoidal function here to make this metric not-so-discontinuous
        case 'max-ratio'
            post_ratio = flow_post ./ flow_max;
            metric(i) = max(post_ratio);
        otherwise
            error('Unknown method');
    end
    if rem(i,100000)==0
        fprintf('Completed %d of %d metric calcs\n',i,n_outages);
    end
end

%keyboard