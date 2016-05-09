function [approved, ps_agent_f] = check_capacity(delta_Pd,delta_Pg,ps_agent_f,ps_agents_local,verbose)
% check if the agents have enough capacity in load/generation to implement
% a solution. In other words, reserve the capacity if it exists
EPS = 1e-4;
C = psconstants;
all_gen_approved = true;
all_sh_approved = true;
for i = 1:length(ps_agents_local)
    % get some information
    ps_agent_loc = ps_agents_local(i);
    bus_id = ps_agent_loc.bus_id;
    ig = (ps_agent_loc.gen(:,C.ge.bus) == bus_id); % find the generation indices
    ish = (ps_agent_loc.shunt(:,C.sh.bus) == bus_id); % find the load indices
    % now find available capacity
    Pg_cap = ps_agent_loc.capacity.Pg;
    sf_cap = ps_agent_loc.capacity.sf;
    Pd_cap = ps_agent_loc.shunt(ish,C.sh.P).*sf_cap;
    if any(Pd_cap < -EPS)
        if all(ps_agent_loc.shunt(ish,C.sh.P) > EPS)
            error('Something is wrong in setting the load capacity.')
        end
        Pd_cap = Pd_cap .* (Pd_cap>EPS);
    end        
    this_ramp_cap = ps_agent_loc.capacity.this_ramp;
    % check if after generation will be positive
    Pg_after = Pg_cap + delta_Pg(ig);
    if all(Pg_after > -EPS)
        gen_approved = true;
    else
        gen_approved = false;
        ps_agent_f.gen(ig,C.ge.P) = Pg_cap; % update agent on the from side
    end
    % check if ramping is allowed
    if -delta_Pg(ig) > this_ramp_cap + EPS
        gen_approved = false;
        ps_agent_f.gen(ig,C.ge.ramp_rate_down) = this_ramp_cap;
    end
    % check if the solving agent knows what the local agent will do
    if max(abs(ps_agent_f.delta_Pg0_all(ig) - ps_agent_loc.delta_Pg)) > EPS
        gen_approved = false;
        ps_agent_f.delta_Pg0_all(ig) = ps_agent_loc.delta_Pg;
    end
    all_gen_approved = all_gen_approved && gen_approved;
    
    ish = (ps_agent_loc.shunt(:,C.sh.bus) == bus_id);
    Pd_after = Pd_cap + delta_Pd(ish);
    if all(Pd_after > -EPS) % same check for loads
        sh_approved = true;
    else
        sh_approved = false;
        ps_agent_f.shunt(ish,C.sh.factor) = sf_cap;
    end
    % check if the solving agent knows what the local agent will do
    if max(abs(ps_agent_f.delta_Pd0_all(ish) - ps_agent_loc.delta_sf .* ps_agent_loc.shunt(ish,C.sh.P))) > EPS
        sh_approved = false;
        ps_agent_f.delta_Pd0_all(ish) = ps_agent_loc.delta_sf .* ps_agent_loc.shunt(ish,C.sh.P);
    end
    all_sh_approved = all_sh_approved && sh_approved;
end
approved = all_gen_approved && all_sh_approved;

% print something
if verbose
    if all_gen_approved
        fprintf('   Generation shedding approved by all agents.\n');
    else
        fprintf(['   Generation shedding not approved by all agents. ', ...
                'Updating agent on bus %d ...\n'], ps_agent_f.bus_id);
    end
    if all_sh_approved
        fprintf('   Load shedding approved by all agents.\n');
    else
        fprintf(['   Load shedding not approved by all agents. ', ...
                'Updating agent on bus %d ...\n'], ps_agent_f.bus_id);
    end
end


