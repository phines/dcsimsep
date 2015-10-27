function ps_agents = set_up_agents(ps,opt)
% initialize agents' models: make a copy from ps, initialize all useful
% information to NaN's and pass it on to the agents
C = psconstants;
ng = size(ps.gen,1);
nd = size(ps.shunt,1);
nbus = size(ps.bus,1);
F = ps.bus_i(ps.branch(:,C.br.from));
T = ps.bus_i(ps.branch(:,C.br.to));

ps_temp = ps;
ps_temp.branch(:,[C.br.Pf,C.br.Pt,C.br.Qf,C.br.Qt]) = NaN;
ps_temp.branch(:,C.br.status) = 1; % the default is that all branches are clsoed
ps_temp.gen(:,[C.ge.Pg,C.ge.Qg]) = NaN;
ps_temp.gen(:,C.ge.status) = 1;
% ps_temp.gen(:,C.ge.ramp_rate_down) = ramp_limits; % don't need this here
ps_temp.shunt(:,C.sh.status) = NaN;
ps_temp.delta_Pg0_all = zeros(ng,1);
ps_temp.delta_Pd0_all = zeros(nd,1);
ps_temp.capacity = [];
ps_temp.delta_Pg = [];
ps_temp.delta_sf = [];
ps_temp.messages = zeros(100,1);
ps_agents(1:nbus) = ps_temp;
adj_mat = sparse([F;T],[T;F],1,nbus,nbus); % the adjacency matrix of the network
% adj_mat^n gives the number of paths with n hops from node i to node j.
% So, the nonzero elements mean that there exists an n-hop path from node
% i to j
A_loc = speye(nbus);
%A_ext = speye(nbus);
%{ 
%extended is the same as local for now
for i = 1:opt.sim.nHopExt
    if i <= opt.sim.nHopLoc
        A_loc = A_loc + adj_mat^i;
    end
    A_ext = A_ext + adj_mat^i;
end
%}
for i = 1:opt.sim.nHopLoc
    A_loc = A_loc + adj_mat^i;
end
% give the local and extended neighborhood bus_i's to each agent
for i = 1:nbus
    ps_agents(i).bus_id = ps.bus(i,C.bu.id); % save bus IDs
    ps_agents(i).loc_nei = find(A_loc(i,:) > 0); % save local neighborhood (bus_i)
    %ps_agents(i).ext_nei = find(A_ext(i,:) > 0); % save extended neighborhood (bus_i)
end



