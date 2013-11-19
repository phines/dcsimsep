function plot_links(~,links,locs,flow,flow_limit,link_status)

base_width = 50;
GREY = [1 1 1]*.5;

m = size(links,1);

for i=1:m
    link = links(i,1:2);
    X = locs(link,1);
    Y = locs(link,2);
    width = max(abs(flow(i))/base_width,0.5);
    if ~link_status(i)
        line( X, Y, 'color', 'k', 'linewidth',5 );
    elseif flow(i)>flow_limit(i)
        line( X, Y, 'color', 'r', 'linewidth', width);
    else
        line( X, Y, 'color', 'g', 'linewidth', width);
    end
end
