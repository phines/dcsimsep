function neighbors = find_neighbors(A,i,distance)
% Find neighbors no more than "distance" hops from "i" in the network
% with adjacency matrix A

n = size(A,1);
is_neighbor = false(n,1);
is_neighbor(i) = true;

for d = 1:distance
    [Ai,~] = find(A(:,is_neighbor));
    is_neighbor(Ai) = true;
end

is_neighbor(i) = false;
neighbors = find(is_neighbor);
