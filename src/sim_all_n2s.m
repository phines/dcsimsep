function sim_all_n2s(set_no,n_subsets,casename,loadprc)
% Simulate all possible n-2 sequences (defaults to the Polish case)
%  this file supports compiled code on the VACC

if ischar(set_no), set_no = str2num(set_no); end
if ischar(n_subsets), n_subsets = str2num(n_subsets); end
if ischar(loadprc), loadprc = str2num(loadprc); end

%% check the inputs and load the data
disp('loading the data');
if strcmp(casename,'polish')
    ps_casename = sprintf('ps_polish_%d',loadprc);
    vars = load('ps_polish_all.mat',ps_casename); % the base case data
    ps = vars.(ps_casename);
end

m = size(ps.branch,1);
all_pairs = list_of_all(m);
np = size(all_pairs,1);

%np = 10; %DEBUG!!!

%% check/process the inputs
if isempty(n_subsets)
    set_no=1;
    subset = (1:np);
else
    % for compiled code...
    if ischar(set_no)
        set_no = str2num(set_no); %#ok<*ST2NM>
        n_subsets = str2num(n_subsets);
    end
    n_per_subset = ceil(np/n_subsets);
    first = (set_no-1)*n_per_subset + 1;
    last  = min(set_no*n_per_subset,np);
    %last = first+5; % DEBUG
    subset = first:last;
end

%% get constants that help us to find the data
%C = psconstants; % tells me where to find my data
% and make sure that ps is updated:
ps = updateps(ps);

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation

%% simulate each n-2 outage
% open a file for the aggregate results
k = 2;
mkdir(sprintf('../%s_results/k2',casename));
fname = sprintf('../%s_results/k%d/bo_sizes_loadprc_%d_%d.csv', ...
    casename,k,loadprc,set_no);
fk = fopen(fname,'w');
if fk<=0, error('could not open file'); end
fprintf(fk,'outages, BO size (MW), BO size (branches)\n');

% open the outage sequence data file
fname = sprintf('../%s_results/k%d/all_k2_loadprc_%d_%d.csv',...
    casename,k,loadprc,set_no);
f = fopen(fname,'w');
if f<=0, error('could not open file'); end

% loop through all of the events
for p = subset
    % simulate the sequence
    fprintf('Simulating n-2 outage %d of %d\n',p,np);
    br_outages = all_pairs(p,:);
    [~,relay_outages,MW_lost] = dcsimsep(ps,br_outages,[],opt);
    fprintf('BO size = %.2f MW\n',MW_lost);
    % save bo size to the aggregate file
    fprintf(fk,'%d,',br_outages);
    fprintf(fk,'%g,%g\n',MW_lost,size(relay_outages,1));

    %fprintf(f,'type (exogenous=0, endogenous=1, stop=-1), time (sec), branch number # (%.4f MW lost)\n', MW_lost);
    % print the exogenous outages
    t = 1.0;
    fprintf(f,'0,%.1f,%d\n',t,br_outages(1));
    fprintf(f,'0,%.1f,%d\n',t,br_outages(2));
    if fk<=0, error('could not open file'); end
    
    % print the endogenous outages
    for i = 1:size(relay_outages,1)
        t = relay_outages(i,1);
        fprintf(f,'1,%.4f,%d\n',t,relay_outages(i,2));
    end
    % print the stop event
    fprintf(f,'-1,%.4f,-1\n',t+100);
end
% close the files
fclose(f);
fclose(fk);
