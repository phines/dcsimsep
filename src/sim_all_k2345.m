% sim_all_k2345
clear all;
addpath('../polish_results/');

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
load k2345results;   % the list of disturbances

m = size(ps.branch,1);

%% simulate each n-k outage
for k = 2:5
    fprintf('Simulating n-%d outages\n',k);
    if     k==2, R = R2;
    elseif k==3, R = R3;
    elseif k==4, R = R4;
    elseif k==5, R = R5;
    end
    % open a file for the aggregate results
    mkdir polish_results
    fname = sprintf('polish_results/results_n_%d.csv',k);
    fk = fopen(fname,'w');
    if fk<=0, error('could not open file'); end
    fprintf(fk,'outages, BO size (MW), BO size (branches)\n');
    
    n_sims = size(R,1);
    for o = 1:n_sims
        % simulate the sequence
        fprintf('Simulating n-%d outage %d of %d\n',k,o,n_sims);
        br_outages = find(R(o,:));
        [isbo,relay_outages,MW_lost] = dcsimsep(ps,br_outages,[],opt);
        fprintf('BO size = %.2f MW\n',MW_lost);
        % save to the aggregate file
        fprintf(fk,'%d,',br_outages);
        fprintf(fk,'%g,%g\n',MW_lost,size(relay_outages,1));
        
        % save the sequence of events to a file
        fname = sprintf('polish_results/k%d_polish_results_%05d.csv',k,o);
        f = fopen(fname,'w');
        fprintf(f,'type (exogenous=0, endogenous=1), time (sec), branch number # (%.4f MW lost)\n', MW_lost);
        % print the exogenous outages
        fprintf(f,'0,1.0,%d\n',br_outages);
        if fk<=0, error('could not open file'); end

        % print the endogenous outages
        for i = 1:size(relay_outages,1)
            fprintf(f,'1,%.4f,%d\n',relay_outages(i,1),relay_outages(i,2));
        end
        fclose(f);
    end
    fclose(fk);
end