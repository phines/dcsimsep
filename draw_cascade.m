function draw_cascade(ps,exo_branches,exo_buses,endo_branches,do_sequence)
% draw the sequence of events for a cascading failure
% exo_branches is exogenous branch outage set
% exo_buses is exogenous bus outage set
% endo_branches is the endogenous outage set, in stages
% if do_sequence is true, then we will plot the whole sequence of events, movie style
% othewise, just a static graph

figure(gcf)
clf

% inputs
if nargin<5
    do_sequence=false;
end

% extract data
C = psconstants;
nodes = ps.bus_i(ps.bus(:,1)); % the first column is always the bus/node number
links = ps.bus_i(ps.branch(:,1:2)); % the first two columns in branch are the end points of the branch
locs  = ps.bus(:,C.bu.locs); % locations
x = locs(:,1);
y = locs(:,2);
n = length(nodes);
m = size(links,1);
link_status = ps.branch(:,C.br.status)~=0;

%constants
GREY = [1 1 1]*.6;

% draw all of the branches
for i=1:m
    link = links(i,1:2);
    X = locs(link,1);
    Y = locs(link,2);
    line( X, Y, 'color', 'r' );
end
axis off; hold on;
% make a title
title([num2str(length(exo_buses)),' buses out:  [',num2str(exo_buses),'];  ',num2str(length(exo_branches)),' branches out:  [',num2str(exo_branches),']']);

% plot the exogenous branch failures
if nargin>1 && ~isempty(exo_branches)
    for br=1:length(exo_branches)
        f = links(exo_branches(br),1); % from end bus
        t = links(exo_branches(br),2); % to end bus
        plot([x(f) x(t)],[y(f) y(t)],'k-','linewidth',4); hold on;
    end
end
link_status(exo_branches) = false;

% only do the following if we are not making a movie.
if ~do_sequence
    % plot the endogenous branch failures
    for br=1:length(endo_branches)
        f = links(endo_branches(br,2),1); % from end bus
        t = links(endo_branches(br,2),2); % to end bus
        plot([x(f) x(t)],[y(f) y(t)],'-','linewidth',4,'color',GREY); hold on;
    end
    link_status(endo_branches(:,2)) = false;
end

% plot the exogenous bus failures
if nargin>2 && ~isempty(exo_buses)
    plot(x(exo_buses),y(exo_buses),'kp','markersize',15,'markerfacecolor','k','markeredgecolor','k');
end

% plot the nodes in groups
plot_nodes(nodes,links(link_status,:),locs);

% put some labels on the graphs
if ~do_sequence
    for br=1:length(exo_branches)
        f = links(exo_branches(br),1); % from end bus
        t = links(exo_branches(br),2); % to end bus
        x_mid = (x(f)+x(t))/2;
        y_mid = (y(f)+y(t))/2;
        text(x_mid,y_mid,'0','HorizontalAlignment','center','BackgroundColor','w');
    end
    
    for br=1:length(endo_branches)
        f = links(endo_branches(br,2),1); % from end bus
        t = links(endo_branches(br,2),2); % to end bus
        x_mid = (x(f)+x(t))/2;
        y_mid = (y(f)+y(t))/2;
        text(x_mid,y_mid,num2str(br),'HorizontalAlignment','center','BackgroundColor','w');
    end
end


% plot the sequence of dependant events
if nargin>3 && do_sequence
    %set(gca,'NextPlot','replacechildren');
    for ti = 1:size(endo_branches,1)
        t = endo_branches(ti,1);
        br_no = endo_branches(ti,2);
        % draw the link outage
        link = links(br_no,1:2);
        X = locs(link,1);
        Y = locs(link,2);
        plot(X,Y,'g-','linewidth',4);
        
        % figure out what the subgraphs are at this point in the simulation
        link_status(br_no) = false;
        % plot the nodes in groups
        plot_nodes(nodes,links(link_status,:),locs);
        [graphNos,nSubGraphs] = find_subgraphs(nodes,links(link_status,:));
        if nSubGraphs>1
            colors = hsv(nSubGraphs);
        else
            colors = 'r';
        end
        for g = 1:nSubGraphs
            % draw the nodes
            color = colors(g,:);
            loc_set = locs(graphNos==g,:);
            plot(loc_set(:,1),loc_set(:,2),'.','MarkerFaceColor',color,'MarkerEdgeColor',color);
        end
        xlabel(['t=',num2str(t)]);
        drawnow
        % record movie frames
        % print something
        fprintf(' t = %6.2f. Branch %4d out. %d islands\n',t,br_no,nSubGraphs);
        %pause
    end
end

return
%%%% end of the fucntion %%%%

% subfunction to just plot the nodes
function plot_nodes(nodes,links,locs)

[graphNos,nSubGraphs] = find_subgraphs(nodes,links);
colors = {'r','b','g','y','w','c','m'};
nc = length(colors);

for g = 1:nSubGraphs
    % draw the nodes
    ci = mod(g,nc);
    if ci==0,ci=length(colors);end
    color = colors{ci};
    loc_set = locs(graphNos==g,:);
    plot(loc_set(:,1),loc_set(:,2),'.',...
        'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',15);
end

    