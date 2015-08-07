load('ps_polish_all','ps_polish_100');
ps = ps_polish_100;
C = psconstants; % tells me where to find my data

nodes = ps.bus(:,1);
links = ps.branch(:,1:2);
A = adjacency(nodes,links);

neighbors1 = find_neighbors(A,17,1)
neighbors2 = find_neighbors(A,14,1)
