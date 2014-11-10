%% constants
GREY = [1 1 1]*.6;
C = psconstants; % tells me where to find my data

%% simulate the cascade
% get constants that help us to find the data
C = psconstants;
% Prepare and run the simulation for the Polish grid
disp('loading the data');
%ps = case300_001_ps;
load case2383_mod_ps;
load BOpairs;

ps = updateps(ps);
m = size(ps.branch,1);

% choose some branch outages
outage_no = 1;
exo_branches = BOpairs(outage_no,:);
exo_buses = [];

% Set some options
opt = psoptions;
opt.verbose = true; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation
opt.sim.fast_ramp_mins = 1;

% run the simulator
disp('Running the simulation and creating the pictures');
dcsimsep2pictures(ps,exo_branches,exo_buses,opt);

return

%% now draw the cascade
% prep the figure
figure(1); clf;
% prep the movie
mov_obj = VideoWriter(sprintf('cascade_%03d',outage_no));
set(mov_obj,'FrameRate',2);
open(mov_obj);
% other stuff
C = psconstants;
nodes = ps.bus_i(ps.bus(:,1)); % the first column is always the bus/node number
links = ps.bus_i(ps.branch(:,1:2)); % the first two columns in branch are the end points of the branch
locs  = ps.bus(:,C.bu.locs); % locations
x = locs(:,1);
y = locs(:,2);
n = length(nodes);
link_status = ps.branch(:,C.br.status)~=0;
flow_limit = ps.branch(:,C.br.rateB);

% draw all of the branches
plot_links(nodes,links,locs,flows(:,1),flow_limit,link_status);
axis off; hold on; axis tight;

% plot the nodes in groups
plot_nodes(nodes,links(link_status,:),locs);

% plot the title and pause
t = 0;
title(sprintf('t = %.2f',t));
%writeVideo(mov_obj,getframe);

% label the exogenous failures
for br=1:length(exo_branches)
    f = links(exo_branches(br),1); % from end bus
    t = links(exo_branches(br),2); % to end bus
    plot([x(f) x(t)],[y(f) y(t)],'k-','linewidth',10); hold on;
end
link_status(exo_branches) = false;
t = 1;
title(sprintf('t = %.2f',t));
writeVideo(mov_obj,getframe);

% plot the sequence of dependant events
%set(gca,'NextPlot','replacechildren');
for ti = 1:size(endo_branches,1)
    cla
    br_no = endo_branches(ti,2);
    % draw the links
    plot_links(nodes,links,locs,flows(:,ti+1),flow_limit,link_status);
    
    % figure out what the subgraphs are at this point in the simulation
    link_status(br_no) = false;
    % plot the nodes in groups
    n_subs = plot_nodes(nodes,links(link_status,:),locs);
    % draw the link outage
    f = links(br_no,1); % from end bus
    t = links(br_no,2); % to end bus
    plot([x(f) x(t)],[y(f) y(t)],'k-','linewidth',10); hold on;

    % print something
    time = endo_branches(ti,1);
    title(sprintf(' t = %6.2f',time));
    fprintf(' t = %6.2f. Branch %4d out. %d islands\n',time,br_no,n_subs);
    if ti<size(endo_branches,1)
        dt = endo_branches(ti+1,2) - t;
    else
        dt = 1;
    end
    drawnow
    % record the movie frames
    %for i = 1:dt
    writeVideo(mov_obj,getframe);
    %end
end


close(mov_obj);

return
%%%% end of the fucntion %%%%

% subfunction to just plot the nodes
