function [graphNos,nSubGraphs] = find_subgraphs(nodes_A,links)
% This function identifies connected components in an undirected graph
%
% usage: [graphNos,nSubGraphs] = find_subgraphs(A)
%  where A is an n x n symetrical incidence matrix 
%or
% usage: [graphNos,nSubGraphs] = find_subgraphs(nodes,links)
%  where nodes is  an n x 1 list of node numbers and links is
%  m x (2+) list of edges (from, to)
%
% The return value is a n x 1 vector of sub-graph numbers, and the
% number of sub-graphs found

[n,dum] = size(nodes_A);
if n==dum
    A = sparse(nodes_A);
else
    nodes = nodes_A;
    n = length(nodes);
    % figure out which links are internal
    internal = ismember(links(:,1),nodes) & ismember(links(:,2),nodes);
    % renumber to get sequential numbering
    e2i = sparse(nodes,1,(1:n)',max(nodes),1);
    % form the incidence matrix
    F = e2i( links(internal,1) );
    T = e2i( links(internal,2) );
    A = sparse([F;T],[T;F],1,n,n) + speye(n);
    %m = size(F,1);
end

grNo = 1;
graphNos = zeros(n,1);
next = 1;
while ~isempty(next)
  included = false(n,1);
  included(next) = true;
  oldLen = 0;
  while sum(included) ~= oldLen
    oldLen = sum(included);
    [Ai,~] = find(A(:,included));
    included(Ai) = true;
  end
  graphNos(included) = grNo;
  grNo = grNo+1;
  next = find(graphNos==0,1);
end

if nargout > 1
  nSubGraphs = max(graphNos);
end

return


