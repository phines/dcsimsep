function drawCascade(branchesout, busesout,brancheslost,gridNoTrace)
figure(gcf)
clf
ps = case300_001_ps;
load case300_locs.txt
nodes = ps.bus(:,1); % the first column is always the bus/node number
links = ps.branch(:,1:2); % the first two columns in branch are the end points of the branch
drawGraph(nodes,links,case300_locs); 
title([num2str(length(busesout)),' buses out [',num2str(busesout),'];  ',num2str(length(branchesout)),' branches out [',num2str(branchesout),']']);
axis on
set(gca,'xtick',[],'ytick',[]);
set(gcf,'color','w')

hold on

plot(case300_locs(busesout,1),case300_locs(busesout,2),'kp','markersize',15,'markerfacecolor','k','markeredgecolor','k');


for br=1:length(branchesout)
    plot(case300_locs(ps.branch(branchesout(br),1:2),1),case300_locs(ps.branch(branchesout(br),1:2),2),'k-','linewidth',4);
end

starttime=find(brancheslost(:,1)>1);
stoptime=find(isnan(brancheslost(:,1)));
for t=starttime:stoptime-1
    plot(case300_locs(ps.branch(brancheslost(t,2),1:2),1),case300_locs(ps.branch(brancheslost(t,2),1:2),2),'g-');
    scatter(case300_locs(nodes,1),case300_locs(nodes,2),10,gridNoTrace(:,t),'filled');
    xlabel(['t=',num2str(brancheslost(t,1))]);
    drawnow
    pause(.1)
end



    