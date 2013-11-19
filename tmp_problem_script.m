ss=[47 63 88 51 89 6 73 9 67 49 113 50 81];
load('ps_50.mat');
ps = updateps(ps);
opts = psoptions;
[is_blackout,relay_outages,MW_lost,p_out,busesout] = dcsimsep(ps,ss,[],opts);