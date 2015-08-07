function [is_sep,sub_grids,n_sub,p_out,busessep] = check_separation(ps,threshold,verbose)%MJE
% returns "NO_SEP" if there is no separation.
% returns "SMALL_SEP" if the fraction of buses in the largest component is less than threshold
% returns "BIG_SEP" if there is major separation

% inputs
if nargin<2, threshold=0.9; end
if nargin<3, verbose=1; end

% extract some data
C = psconstants;
n = size(ps.bus,1);

% set up the outputs
NO_SEP = 0;
BIG_SEP = 2;
SMALL_SEP = 1;

% find the subgraphs
br_status = (ps.branch(:,C.br.status)~=0);
[sub_grids,n_sub] = find_subgraphs(ps.bus(:,1),ps.branch(br_status,1:2));

% if there are no subgraphs, 
if n_sub == 1
   is_sep = NO_SEP;
   p_out=0; %MJE
   busessep=[]; %MJE
   return
else
   is_sep = SMALL_SEP;
   if verbose
       fprintf('System is divided into %d subgrids\n',n_sub);
   end
end

% check for major separation
sub_grid_size = zeros(n_sub,1);
for g = 1:n_sub
   subset = (sub_grids==g);
   sub_grid_size(g) = sum(subset);
end

% find the largest component
[s,maxsubset] = max(sub_grid_size); %MJE
busessep=find(sub_grids~=maxsubset); %find which buses were separated %MJE

if s/n<threshold
   is_sep = BIG_SEP;
   if verbose
       fprintf('Major separation found\n');
   end
end
p_out=1-s/n; %return proportion of the buses that are separated %MJE
