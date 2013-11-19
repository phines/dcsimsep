% test_dcpf

% load data
C = psconstants;
ps = case6_ps;
% run the base case
ps = dcpf(ps);
printps(ps);
disp('Above are the base case results');
pause

% separate out bus 1
ps = case6_ps;
ps.branch(:,C.br.status) = 1;
ps.branch(1:3,C.br.status) = 0;
ps = dcpf(ps);
printps(ps);
disp('Above are the results after separating out bus 1');
pause

% separate out buses 1 and 2
ps = case6_ps;
ps.branch(:,C.br.status) = 1;
ps.branch([1 3 5 10],C.br.status) = 0;
ps = dcpf(ps);
printps(ps);
disp('Above are the results after separating out buses 1 and 2');

% separate out buses 4 and 5
ps = case6_ps;
ps.branch(:,C.br.status) = 1;
ps.branch([2 3 5 6 8 11],C.br.status) = 0;
ps = dcpf(ps);
printps(ps);
disp('Above are the results after separating out buses 4 and 5');
