function PTDF = get_ptdf(ps,use_matpower)
% gets a power transfer distribution factor matrix
C = psconstants;
if nargin<2
    use_matpower = false;
end

if use_matpower
    ps_i = ps; % make ps have consecutive bus numbers
    ps_i.bus(:,C.bu.id) = full(ps.bus_i(ps.bus(:,C.bu.id)));
    ps_i.branch(:,C.br.from) = full(ps.bus_i(ps.branch(:,C.br.from)));
    ps_i.branch(:,C.br.to) = full(ps.bus_i(ps.branch(:,C.br.to)));
    
%     ps_i.branch(:,C.br.B) = 0;
%     ps_i.branch(:,C.br.tap) = 1;
%     ps_i.branch(:,C.br.shift) = 0;
    
    PTDF = makePTDF(ps_i.baseMVA, ps_i.bus, ps_i.branch);
    %LODF = makeLODF(ps_i.branch, PTDF);
else
    % get some data
    n = size(ps.bus,1);
    m = size(ps.branch,1);
    br_st = ps.branch(:,C.br.status);
    F = full(ps.bus_i(ps.branch(:,1)));
    T = full(ps.bus_i(ps.branch(:,2)));
    inv_X = br_st./ps.branch(:,C.br.X);

    % build the B bus matrix
    B = sparse(F,T,-inv_X,n,n) + ...
        sparse(T,F,-inv_X,n,n) + ...
        sparse(T,T,+inv_X,n,n) + ...
        sparse(F,F,+inv_X,n,n);
    
    % make the PTDF
    PTDF = zeros(m,n);
    Bf   = sparse(1:m,F,+inv_X,m,n) + ...
           sparse(1:m,T,-inv_X,m,n);
        
    nonref = true(n,1);
    ref_bus = find(ps.bus(:,C.bu.type)==3);
    if length(ref_bus) == 1
        nonref(ref_bus) = false;
        PTDF(:,nonref) = full(Bf(:, nonref) / B(nonref, nonref));
    else
        error('More than 1 refrence bus detected!')
    end
% The following is my implementation of the 2-page paper on direct 
% calculation of LODf, which turns out to be same as above!

%     LODF = zeros(m,m);
%     for i=1:m
%         sai = sparse(F(i),1,1,n,1) + ...
%               sparse(T(i),1,-1,n,1);
%         PTDFmo = PTDF(:,nonref)*sai(nonref);
%         h = B(nonref,nonref)\sai(nonref);
%         PTDFoo = inv_X(i)*sai(nonref)'*h;  
%         LODF(:,i) = PTDFmo/(1-PTDFoo);
%     end
%     LODF = LODF - diag(diag(LODF)) - eye(m, m);
end
