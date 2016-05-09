function ps = implement_solution(ps,ps_agents,verbose)
C = psconstants;
EPS = 1e-4;
ng = size(ps.gen,1); nd = size(ps.shunt,1);
all_dPg = zeros(ng,1); all_dPd = zeros(nd,1); % keep track to print output
for i = 1:length(ps_agents)
    % find gen indices and implement delta_Pg
    this_agent = ps_agents(i);
    bus_id = this_agent.bus_id;
    % find the generation indices and update ps.gen 
    ig = (this_agent.gen(:,C.ge.bus) == bus_id);
    ps.gen(ig,C.ge.P) = ps.gen(ig,C.ge.P) + this_agent.delta_Pg;
    all_dPg(ig) = this_agent.delta_Pg;
    % find the shunt indices and update ps.shunt
    ish = (this_agent.shunt(:,C.sh.bus) == bus_id);
    ps.shunt(ish,C.sh.factor) = ps.shunt(ish,C.sh.factor) + this_agent.delta_sf;
    all_dPd(ish) = this_agent.delta_sf .* ps.shunt(ish,C.sh.P);
end
idx = abs(ps.gen(:,C.ge.P)) < EPS;
ps.gen(idx,C.ge.P) = 0;
idx = abs(ps.shunt(:,C.sh.factor)) < EPS;
ps.shunt(idx,C.sh.factor) = 0;
if verbose
    print_load_shedding_results(ps,all_dPg,all_dPd)
end


