function [nodes,links,node_locs] = build_reg_grid(n,m)
% Builds a regular grid, with nodes/links randomly removed
% until the size of the graph is n-nodes and m-links
% Outputs:
%  nodes - a list of nodse
%  links - a list of edge end points
%  node_locs - x/y coordinates for the nodes (for plotting)

% get the lenght of a side
len = ceil(sqrt(n));

x = []; y = [];
f = []; t = [];
for i = 1:len
    % new node IDs
    id = (1:len)' + (i-1)*len;
    % node locations
    x_ = ones(len,1)*i;
    y_ = (1:len)';
    x = [x;x_]; %#ok<*AGROW>
    y = [y;y_];
    % links from left to right
    f1 = id(1:len-1);
    t1 = id(2:len);
    % links from top to bottom
    if i<len
        f2 = id(1:len);
        t2 = id(1:len) + len;
    end
    f = [f;f1;f2];
    t = [t;t1;t2];
end

% Build the info about the graph
nodes = union(f,t);
links = [f t];
node_locs = [x y];

% Remove nodes randomly until we have n nodes
n_hat = length(nodes);
while n_hat>n
    % choose a random number
    r = randi([1 n_hat]);
    % delete the node
    node_name = nodes(r);
    nodes_tmp = nodes;
    nodes_tmp(r) = [];
    % delete links that include that node
    links_tmp = links;
    dead_links = ismember(links(:,1),node_name) | ismember(links(:,2),node_name);
    if any(dead_links)
        links_tmp(dead_links,:) = [];
    end
    % test for connectedness
    [~,n_sub] = find_subgraphs(nodes_tmp,links_tmp);
    if n_sub==1
        nodes = nodes_tmp;
        links = links_tmp;
        node_locs(r,:) = [];
        n_hat = length(nodes);
    end
end

% Remove links randomly until the number of links is m
m_hat = size(links,1);
if m_hat<m
    error('Not possible to build a regular graph of this size');
end
% randomly remove links
while m_hat > m
    % choose a random number
    r = randi([1 m_hat]);
    % delete the link
    links_tmp = links;
    links_tmp(r,:) = [];
    % test for connectedness
    [~,n_sub] = find_subgraphs(nodes,links_tmp);
    if n_sub==1
        links = links_tmp;
        m_hat = size(links,1);
    end
end
