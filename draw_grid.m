clear;clc;
n = 100; % n x n grid

x = []; y = [];
f = []; t = [];
for i = 1:n
    % new node IDs
    id = (1:n)' + (i-1)*n;
    % node locations
    x_ = ones(n,1)*i;
    y_ = (1:n)';
    x = [x;x_];
    y = [y;y_];
    % links from left to right
    f1 = id(1:n-1);
    t1 = id(2:n);
    % links from top to bottom
    if i<n
        f2 = id(1:n);
        t2 = id(1:n) + n;
    end
    f = [f;f1;f2];
    t = [t;t1;t2];
end

f(500:1000) = [];
t(500:1000) = [];

% fnDrawNetwork(f,t,x,y);
G = sparse(f,t,ones(1,length(f)),max(t),max(t));
h = view(biograph(G,[],'ShowWeights','on'));
[S, C] = graphconncomp(G);   % find connected components
colors = jet(S);
for i = 1:numel(h.Nodes)
h.Nodes(i).Color = colors(C(i),:);
end

