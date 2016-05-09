clear all;

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation

%% Load the data
disp('loading the data');
load case2383_mod_ps; % the list of 
load pairdata4paul;   % the list of disturbances

m = size(ps.branch,1);

%% simulate each branch outage
n_sims = size(BOpairs,1);
for o = 1:n_sims
    fprintf('Simulating outage %d of %d\n',o,n_sims);
    br_outage = BOpairs(o,:);
    [isbo(o),relay_outage{o},MW_lost(o)] = dcsimsep(ps,br_outage,[],opt); %#ok<SAGROW>
    fprintf('BO size = %.2f MW\n',MW_lost(o)');
end
save n_minus_2_results;


%% save the results in a file
cd polish_results;
load n_minus_2_results;

for o = 1:n_sims
    % type, time, component number
    fname = sprintf('polish_results_%03d.csv',o)
    f = fopen(fname,'w');
    fprintf(f,'type (exogenous=0, endogenous=1), time (sec), branch number # %.4f MW lost\n', MW_lost(o));
    % print the exogenous outages
    fprintf(f,'0,1.0,%d\n',BOpairs(o,1));
    fprintf(f,'0,1.0,%d\n',BOpairs(o,2));
    % print the endogenous outages
    for i = 1:size(relay_outage{o},1)
        fprintf(f,'1,%.4f,%d\n',relay_outage{o}(i,1),relay_outage{o}(i,2));
    end
    fclose(f);
end
cd ..
