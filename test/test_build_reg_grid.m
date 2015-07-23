% build a regular graph
clear; clc;
n = 2383;
m = 2886;

[nodes,links,node_locs] = build_reg_grid(n,m);

[~, ~, links_idx] = unique(links);
nodes_f = links_idx(1:m);
nodes_t = links_idx(m+1:2*m);
links = [nodes_f nodes_t];

[~, ~, nodes_idx] = unique(nodes);
nodes = nodes_idx;
 
figure(1); clf;
draw_graph(nodes,links,node_locs);
