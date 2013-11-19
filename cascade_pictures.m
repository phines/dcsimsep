% Show some interesting cascading failure paths

%% load the data
C = psconstants;
load case2383_mod_ps;
load pairdata4paul
pairs = BOpairs;

ps = dcpf(ps);
flow = ps.branch(:,C.br.Pf);

opt = psoptions;
opt.verbose = false;
m = size(ps.branch,1);
do_pause = 0;

%% try a specific pair
%opt.verbose = 1;
opt.sim.stop_threshold = 0;
pair = pairs(1,:);
[is_blackout,relay_outages,MW_lost] = dcsimsep(ps,pair,[],opt);
figure(1); clf;
draw_cascade(ps,pair,[],relay_outages);
savepng(['cascade' num2str(i)],300);

%% run each contingency
%{
np = size(pairs,1);
p_out = zeros(np,1);
is_blackout = zeros(np,1);
relay_outages = cell(np,1);
MW_lost = zeros(np,1);

for i = 1:size(pairs,1)
    pair = pairs(i,:);
    [is_blackout(i),relay_outages{i},MW_lost(i),p_out(i)] = dcsimsep(ps,pair,[],opt);
%    figure(1); clf;
%    draw_cascade(ps,pair,[],relay_outages);
%    savepng(['cascade' num2str(i)],300);
end
%}

%%
%{
for i = 282%i = 1:size(pairs,1)
    pair = pairs(i,:);
    figure(1); clf;
    draw_cascade(ps,pair,[],relay_outages{i});
    savepng(['cascade' num2str(i)],300);
end
%}


%% draw the system with the BO 
branch_set = [pairs(:,1); pairs(:,2)];

branch_id = unique(branch_set);
count = zeros(size(branch_id));
for i = 1:length(un)
    br = branch_id(i);
    count(i) = sum(branch_set==br);
end

%rows = [branch_id,count];
%rows = sortrows(rows);

%% do the plot
figure(2); clf;
GREY = [1 1 1]*.6;
nodes = ps.bus_i(ps.bus(:,1)); % the first column is always the bus/node number
links = ps.bus_i(ps.branch(:,1:2)); % the first two columns in branch are the end points of the branch
locs  = ps.bus(:,C.bu.locs); % locations
x = locs(:,1);
y = locs(:,2);
m = size(links,1);

axis off; hold on;

% draw all of the branches
for i=1:m
    link = links(i,1:2);
    X = locs(link,1);
    Y = locs(link,2);
    wid = max(abs(flow(i))/100,0.1);
    line( X, Y, 'color', 'b', 'linewidth',wid ); hold on;        
end

% plot the nodes
plot(x,y,'.',...
    'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15);

% highlight the top n branches
index = find( count>=10 )';
for i = index
    br = branch_id(i);
    c  = count(i);
    f = links(links(br,2),1); % from end bus
    t = links(links(br,2),2); % to end bus
    wid = max(abs(flow(br))/100,0.1);
    plot([x(f) x(t)],[y(f) y(t)],'y-','linewidth',wid+5); hold on;
    plot([x(f) x(t)],[y(f) y(t)],'B-','linewidth',wid);
end
% do the text
for i = index
    br = branch_id(i);
    c  = count(i);
    f = links(links(br,2),1); % from end bus
    t = links(links(br,2),2); % to end bus
    x_mid = (x(f)+x(t))/2 + 50;
    y_mid = (y(f)+y(t))/2;
    %text(x_mid,y_mid,num2str(c),...
    %    'HorizontalAlignment','left','BackgroundColor','y','FontSize',14);
    if br==96
        text(x_mid,y_mid,'this is a lot of text');
    end
end
