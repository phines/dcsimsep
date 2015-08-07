function A_blkdiag = make_block_diag(A,n,extra_col,ncol)
% make a sparse block diagonal matrix from A being repeated n times on the 
% diagonal
% extra_col: the extra columns with zeros on the left hand side of the 
% final matrix (useful to build optimization matrices), default: 0
% ncol: total number of columns, defaul: ny*n

[r,c,v] = find(A);
[nx,ny] = size(A);
if nargin < 3, extra_col = 0; end
if nargin < 4, ncol = ny*n; end
% set indices for rows
new_r = bsxfun(@plus,r(:),0:nx:(n-1)*nx);
% set indices for columns
new_c = bsxfun(@plus,c(:),0:ny:(n-1)*ny) + extra_col;
% new values
new_v = repmat(v(:),n,1);

A_blkdiag = sparse(new_r,new_c,new_v,nx*n,ncol);
