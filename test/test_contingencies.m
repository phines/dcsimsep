clear all;
randseed(1); % deterministic randomness

%% load the data
load pairdata4paul;
load case2383_mod_ps;
ps = dcpf(ps);

m = size(ps.branch,1);

% my method
%metric_BO_PCP   = get_contingency_metrics(ps,BOpairs,'PCP');
%metric_NOBO_PCP = get_contingency_metrics(ps,NOBOpairs,'PCP');
% wollenberg method
metric_BO_Woll   = get_contingency_metrics(ps,BOpairs,'wollenberg');
metric_NOBO_Woll = get_contingency_metrics(ps,NOBOpairs,'wollenberg');
%metric_rand_Woll = get_contingency_metrics(ps,random_pairs,'wollenberg');

%% compute FPR/FNR
%estimate the false positive and false negative rate at different values of 
% the wollenberg metric (w)
% check to see if the random pairs of BOs
load all_woll_results

n_all = length(metric_all_Woll);
n_bo = sum(is_bo_all);
n_nbo = n_all-n_bo;
dw_list = 10.^(-3:.1:2.7);
FPR = zeros(1,length(dw_list));
FNR = zeros(1,length(dw_list));
for i = 1:length( dw_list )
    dw = dw_list(i);
    is_gt_dw = metric_all_Woll>dw;
    n_gt_dw = sum(is_gt_dw);
    n_lt_dw = n_all-n_gt_dw;
    n_bo_gt_dw  = sum(is_bo_all(is_gt_dw));  % this is the number of true positives
    n_nbo_gt_dw = (n_gt_dw - n_bo_gt_dw);    % this is the number of false positives
    n_bo_lt_dw  = sum(is_bo_all(~is_gt_dw)); % this is the number of false negatives
    n_nbo_lt_dw = n_lt_dw - n_bo_lt_dw;      % this is the number of true negatives
    FPR(i) = n_nbo_gt_dw / n_nbo;
    FNR(i) = n_bo_lt_dw / n_bo;
end


%% Wollenberg figure

figure(1); clf;
%{
subplot(2,1,1);
set(gca,'FontSize',14);
hold on;
plot_ccdf(metric_BO_Woll,  min(metric_BO_Woll),  'k-');
plot_ccdf(metric_NOBO_Woll,min(metric_NOBO_Woll),'k--');
plot_ccdf(metric_all_Woll,min(metric_all_Woll),'k-.');
set(gca,'xscale','log');
set(gca,'yscale','log');
axis tight;
axis([.01 max(metric_all_Woll) 10^-7 1]);
set(gca,'ytick',10.^(-7:0));
h = legend('BO','Non-BO, matched flow','all');
set(h,'EdgeColor',[1 1 1]);
%}

%subplot(2,1,2); cla;
set(gca,'FontSize',14);
hold on;
plot(dw_list,FPR,'k-');
plot(dw_list,FNR,'k--');
set(gca,'yscale','log');
set(gca,'xscale','log');
axis tight
h = legend('False Positive Rate','False Negative Rate');
set(h,'EdgeColor',[1 1 1]);
axis([.01 max(metric_all_Woll) 10^-7 1]);
set(gca,'ytick',10.^(-7:2:0));
set(gca,'yminortick','on');
%xlabel('Change in Performance Index (\Delta PI_{\overline{r,s})');


return

%% PCP figure
figure(2); clf;
subplot(2,1,1)
hold on;
plot_ccdf(metric_BO_PCP,min(metric_BO_PCP),'k-');
plot_ccdf(metric_NOBO_PCP,min(metric_NOBO_PCP),'k--');
title('PCP');
h = legend('Blackouts','non-blackouts');
set(h,'EdgeColor',[1 1 1]);

return
% junk
% histograms
figure(1); clf;
subplot(2,1,1);
hist([metric_BO_PCP metric_NOBO_PCP],40);
legend('BO-PCP','No BO-PCP');
subplot(2,1,2);
hist([metric_BO_Woll metric_NOBO_Woll],40);
h = legend('BO-PCP-Woll','No BO-PCP-Woll');
set(h,'EdgeColor',[1 1 1]);

% BCP figure

%


%subplot(2,1,2);
%hist(metric_NOBO,20);
