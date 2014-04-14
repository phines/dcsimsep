% build a regular graph

n = 2383;
m = 2886;

[nodes,links,node_locs] = build_reg_grid(n,m);
 
figure(1); clf;
draw_graph(nodes,links,node_locs);
