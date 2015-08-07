function nSubGraphs = plot_nodes(nodes,links,locs)


[graphNos,nSubGraphs] = find_subgraphs(nodes,links);
map = colormap('hsv');
randseed(1); % deterministic randomness
R = randperm(size(map,1));

colors = map(R,:);
nc = length(colors);

for g = 1:nSubGraphs
    % draw the nodes
    ci = mod(g,nc);
    if ci==0,ci=length(colors);end
    color = colors(ci,:);
    loc_set = locs(graphNos==g,:);
    plot(loc_set(:,1),loc_set(:,2),'.',...
        'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',25);
end
