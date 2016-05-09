function printps(ps)
%% prints power system data in the ps structure
C = psconstants;
j = 1i;

%% figure out the generation and load
n = size(ps.bus,1);
D = ps.bus_i(ps.shunt(:,1)); % shunt locations
Sd = (ps.shunt(:,C.sh.P) + j*ps.shunt(:,C.sh.Q)).*ps.shunt(:,C.sh.factor);
Sd_bus = sparse(D,1,Sd,n,1);

G = ps.bus_i(ps.gen(:,1)); % generator locations
Sg = (ps.gen(:,C.ge.P) + j*ps.gen(:,C.ge.Q)).*ps.gen(:,C.ge.status);
Sg_bus = sparse(G,1,Sg,n,1);

%% Print bus data
disp('------------------------------ Bus Data -----------------------------');
disp('   ID  |V|       /_V   Type  Pd(MW)   Qd(MVAr)  Pg(MW)   Qg(MVAr)');
disp('---------------------------------------------------------------------');

for i = 1:n
    Sg = full(Sg_bus(i));
    Sd = full(Sd_bus(i));
    
    type = ps.bus(i,C.bu.type);
    if type==C.REF
        type_str = 'Ref';
    elseif type==C.PV
        type_str = 'PV';
    elseif type==C.PQ
        type_str = 'PQ';
    else
        type_str = 'Oth';
    end
    
    fprintf('%5d %7.4f %8.3f %3s',ps.bus(i,[C.bu.id C.bu.Vmag C.bu.Vang]),type_str);
    fprintf('%8.2f+%8.2fj %8.2f+%8.2fj\n',real(Sd),imag(Sd),real(Sg),imag(Sg));
end

%% Print branch data
nbr = size(ps.branch,1);
disp('----------------------------- Branch Data ------------------------------------------------');
disp(' ID  From    To  st   Pf(MW)  Qf(MVAr)   Pt(MW) Qt(MVAr)   Rate A  |If|(A)  |It|(A)  ');
disp('------------------------------------------------------------------------------------------');
for i = 1:nbr
    from = ps.branch(i,1);
    to   = ps.branch(i,2);
    Sft = ps.branch(i,C.br.Pf:C.br.Qt);
    baseV_lg_f = ps.bus(ps.bus_i(from),C.bu.baseKV) * 1000 / sqrt(3);
    baseV_lg_t = ps.bus(ps.bus_i(to  ),C.bu.baseKV) * 1000 / sqrt(3);
    baseIf = ps.baseMVA*1e6 / baseV_lg_f / 3;
    baseIt = ps.baseMVA*1e6 / baseV_lg_t / 3;
    If = ps.branch(i,C.br.Imag_f) * baseIf;
    It = ps.branch(i,C.br.Imag_t) * baseIt;
    rateA = ps.branch(i,C.br.rateA);
    fprintf('%3d %5d %5d %3d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n', ...
        [i from to ps.branch(i,C.br.status) Sft rateA If It]);
end

%% Print gen data
disp('------------------ Generator Data --------------------');
disp('  Bus st   Pmin P(MW) Pmax  Qmin  Q(MW)  Qmax ');
disp('------------------------------------------------------');
fprintf('%5d  %1d %4.0f %7.2f %4.0f %4.0f %7.2f %4.0f\n',ps.gen(:,[1 C.ge.status C.ge.Pmin C.ge.P C.ge.Pmax C.ge.Qmin C.ge.Q C.ge.Qmax ])');
%% Print system data
disp('------------------- System Data ----------------------');
Pg = sum(ps.gen(:,C.ge.P) .* ps.gen(:,C.ge.status));
Qg = sum(ps.gen(:,C.ge.Q) .* ps.gen(:,C.ge.status));
Pd = sum(ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.status));
Qd = sum(ps.shunt(:,C.sh.Q) .* ps.shunt(:,C.sh.status));
fprintf(' Active power losses = %.3f MW, %.2f%%\n',Pg - Pd, (Pg - Pd)/Pg*100 );
fprintf(' Reactive power losses = %.3f MVAr\n',Qg-Qd);
disp('------------------------------------------------------');

 

