% test the mp2ps function

mpc = case3375wp;
ps = mp2ps(mpc);
ps = updateps(ps);
disp('First power flow');
ps = dcpf(ps);
writeps(ps,'case3375wp_ps');
ps = case3375wp_ps;
disp('Second power flow');
ps = dcpf(ps);
pause
printps(ps);
