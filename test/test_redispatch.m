% test_redispatch
clear all;
% load data
C = psconstants;
%{
ps = updateps(case6_ps);
% run the base case
ps = dcpf(ps);
printps(ps);
disp('Above are the base case results');
pause

% separate out bus 1
ps = updateps(case6_ps);
ps.branch(:,C.br.status) = 1;
ps.branch(1:3,C.br.status) = 0;
br_st = ps.branch(:,C.br.status)==1;
[sub_grids,~] = findSubGraphs(ps.bus(:,1),ps.branch(br_st,1:2));
[Pg,ge_status,sh_factor] = redispatch(ps,sub_grids);
ps.gen(:,C.ge.P) = Pg;
ps.gen(:,C.ge.status) = ge_status;
ps.shunt(:,C.sh.factor) = sh_factor;
ps = dcpf(ps,sub_grids);
printps(ps);
disp('Above are the results after separating out bus 1');
pause

% divide the system in half
ps = updateps(case6_ps);
ps.branch(:,C.br.status) = 1;
ps.branch([1 5 6 8 11],C.br.status) = 0;
br_st = ps.branch(:,C.br.status)==1;
[sub_grids,n_sub] = findSubGraphs(ps.bus(:,1),ps.branch(br_st,1:2));
[Pg,ge_status,sh_factor] = redispatch(ps,sub_grids);
ps.gen(:,C.ge.P) = Pg;
ps.gen(:,C.ge.status) = ge_status;
ps.shunt(:,C.sh.factor) = sh_factor;
ps = dcpf(ps,sub_grids);
printps(ps);
disp('Above are the results after separating out buses 1, 4, 5');
pause
%}

% cut the ramp rate down
ps = updateps(case6_ps);
ps.branch(:,C.br.status) = 1;
ps.gen(:,C.ge.Pmax) = 100;
ps.gen(:,C.ge.P) = 100;
br_st = ps.branch(:,C.br.status)==1;
[sub_grids,n_sub] = findSubGraphs(ps.bus(:,1),ps.branch(br_st,1:2));
ramp_rate = ones(3,1);
[Pg,ge_status,sh_factor] = redispatch(ps,sub_grids,ramp_rate);
ps.gen(:,C.ge.P) = Pg;
ps.gen(:,C.ge.status) = ge_status;
ps.shunt(:,C.sh.factor) = sh_factor;
ps = dcpf(ps,sub_grids);
printps(ps);
disp('Above are the results after separating out buses 1, 4, 5, with Pmax=100');

% reduce the generator limits
ps = updateps(case6_ps);
ps.branch(:,C.br.status) = 1;
ps.branch([1 5 6 8 11],C.br.status) = 0;
ps.gen(:,C.ge.Pmax) = 100;
ps.gen(:,C.ge.P) = 100;
br_st = ps.branch(:,C.br.status)==1;
[sub_grids,n_sub] = findSubGraphs(ps.bus(:,1),ps.branch(br_st,1:2));
[Pg,ge_status,sh_factor] = redispatch(ps,sub_grids);
ps.gen(:,C.ge.P) = Pg;
ps.gen(:,C.ge.status) = ge_status;
ps.shunt(:,C.sh.factor) = sh_factor;
ps = dcpf(ps,sub_grids);
printps(ps);
disp('Above are the results after separating out buses 1, 4, 5, with Pmax=100');
