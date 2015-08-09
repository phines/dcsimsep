function [ps, MW_lost, imbalance] = dist_control(ps_agents,ps,sub_grids,ramp_rate,opt)
% distributed emergency control
C = psconstants;
EPS = 1e-4;
nbus = size(ps.bus,1);
ramp_dt = opt.sim.fast_ramp_mins*60;
ramp_limits = ramp_rate * ramp_dt;
verbose = opt.verbose;
Pd0 = ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.factor);
G = ps.bus_i(ps.gen(:,C.ge.bus));
D = ps.bus_i(ps.shunt(:,C.sh.bus));
imbalance = 0;

for i = 1:nbus
    % all agents collect data from the ps structure (measurements)
    collect_data(i,ps)
    % reset capacity and ramp rates for each agent
    reset_agent(i,ramp_limits);
end

% each agent checks for overloads on the lines connected to it on its from side
for bus_i_f = 1:nbus
    % find lines connected to the bus (on the from side)
    [br_ID_f, ~] = find(ps_agents(bus_i_f).branch(:,C.br.f) == ps_agents(bus_i_f).bus_id); % find branches that this bus is on the from side
    flow_max_agent = ps_agents(bus_i_f).branch(br_ID_f,opt.opf.contg_rate);
    idx = (abs(ps_agents(bus_i_f).branch(br_ID_f,C.br.Pf)) > flow_max_agent + EPS);
    if sum(idx) > 0
        br_over_id = br_ID_f(idx);
        bus_i_t = ps_agents(bus_i_f).bus_i(ps_agents(bus_i_f).branch(br_over_id,C.br.to));
        bus_i_t = full(bus_i_t);
        if verbose
            bus_id_f = ps_agents(bus_i_f).bus(bus_i_f,1);
            bus_id_t = ps_agents(bus_i_f).bus(bus_i_t,1);
            fprintf('Bus %d found overload on line(s): ', bus_id_f);
            fprintf('%d(%d-%d) ', ...
                [br_over_id';repmat(bus_id_f,1,length(br_over_id));bus_id_t']);
            fprintf('\n');
        end
        % pulling information from local and extended neighborhoods
        local_nei = unique([ps_agents([bus_i_f;bus_i_t]).loc_nei]);
        local_nei = local_nei(:)'; % make sure local_nei is a row vector
        for i = local_nei
            % prepare the local message
            msg_loc = prepare_msg_loc(ps_agents(i));
            % now send the message to the from bus (in charge of
            % solving)
            agent_send_msg(msg_loc,'local',bus_i_f)
        end
        extend_nei = unique([ps_agents([bus_i_f;bus_i_t]).ext_nei]);
        for i = extend_nei
            % prepare the extended message
            msg_ext = prepare_msg_ext(ps_agents(i));
            % agents in the extended neighborhood send outage data to the from bus
            agent_send_msg(msg_ext,'extended',bus_i_f)
        end
        OptVarBusID = ps_agents(bus_i_f).bus(local_nei,1); % bus IDs of the local neighborhood
        approved = false; % this will change to true if the optimization is solved and there is enough capacity for implementation
        continue_flag = false;
        while ~approved
            if opt.sim.use_mpc
                OptVarBus_I = ps_agents(bus_i_f).bus_i(OptVarBusID);
                [delta_Pd, delta_Pg] = mpc_smp_solve_lp(ps_agents(bus_i_f),OptVarBus_I,ramp_limits,opt);
            else
                [delta_Pd,delta_Pg] = emergency_control_dec(ps_agents(bus_i_f),OptVarBusID,opt);
            end
            % do some checks 
            if abs(sum(delta_Pd)-sum(delta_Pg)) > 1e-3
                error('load shedding is not equal to generation reduction.')
            end
            if sum(delta_Pd)==0 || sum(delta_Pg)==0
                continue_flag = true;
                if verbose
                    fprintf('   zero load shedding/generator reduction found.\n')
                end
                break
            end
            % make sure delta_Pd and delta_Pg are only non-zero on the
            % local neighborhood
            if ~isempty(setdiff(D(delta_Pd ~= 0), local_nei)) || ...
                    ~isempty(setdiff(G(delta_Pg ~= 0), local_nei))
                error('load/generation shedding outside local neighborhood.')
            end
            % check to see if the local agents have capacity
            [approved, ps_agents(bus_i_f)] = check_capacity(delta_Pd,delta_Pg,ps_agents(bus_i_f),ps_agents(local_nei),verbose);
        end
        if continue_flag 
            continue 
        end
        % now implement solution on ps and change capacity for
        % agents
        if verbose  
            fprintf('  sending solution from bus %d to the local neighborhood...\n', ...
                ps_agents(bus_i_f).bus_id); 
        end
        update_agent_data(local_nei,delta_Pd,delta_Pg);
    end
end
% implement all solutions from each agent on ps
ps = implement_solution(ps,ps_agents,verbose);
Pd = ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor);
Pg = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status);
mis = sum(Pg) - sum(Pd);
if mis > EPS
    % record and do a rebalance
    imbalance = imbalance + mis; 
    ps = rebalance(ps,sub_grids,ramp_limits,verbose);
end

% run power flow on ps 
ps = dcpf(ps,sub_grids);
Pg = ps.gen(:,C.ge.Pg) .* ps.gen(:,C.ge.status);
Pg_max = ps.gen(:,C.ge.Pmax) .* ps.gen(:,C.ge.status) + EPS;
Pg_min = ps.gen(:,C.ge.Pmin) .* ps.gen(:,C.ge.status) - EPS;
if any( round(Pg)<round(Pg_min-EPS) | round(Pg)>round(Pg_max+EPS) )
    error('Pg is out of bounds'); 
end

Pd = ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.factor);
if abs(sum(Pd)-sum(Pg)) > EPS
    error('Something is wrong.')
else
    MW_lost = sum(Pd0) - sum(Pd);
end
    
    function agent_send_msg(msg,msg_type,nei)
        % ps_agents_sub is a subset of agents that receive the message. The 
        % extended neighborhood contains the local neighborhood in this
        % implementation, so the message for extended neighborhood goes to the
        % local nodes as well.
        for j = nei
            br_id = msg.branch.status(:,1);
            % update the status of the lines connected to the sending bus (this is
            % the same for both local and extended neighborhood)
            ps_agents(j).branch(br_id,C.br.status) = msg.branch.status(:,2);
            if strcmp(msg_type,'local')
                br_id_f = msg.branch.Pf(:,1);
                br_id_t = msg.branch.Pt(:,1);
                % update flows
                ps_agents(j).branch(br_id_f,C.br.Pf) = msg.branch.Pf(:,2);
                ps_agents(j).branch(br_id_t,C.br.Pt) = msg.branch.Pt(:,2);
                % update the other end flow (this is not part of message passing,
                % but because this is a dc model, giving info about Pf gives info
                % also about Pt and vice versa)
                ps_agents(j).branch(br_id_t,C.br.Pf) = -msg.branch.Pt(:,2);
                ps_agents(j).branch(br_id_f,C.br.Pt) = -msg.branch.Pf(:,2);

                % update loads
                i_sh = msg.shunt.factor(:,1);
                ps_agents(j).shunt(i_sh,C.sh.factor) = msg.shunt.factor(:,2);

                % update generations
                ig = msg.gen.Pg(:,1);
                ps_agents(j).gen(ig,C.ge.Pg) = msg.gen.Pg(:,2);
                ps_agents(j).gen(ig,C.ge.status) = msg.gen.status(:,2);
            elseif strcmp(msg_type,'extended')
                % the status of branches are already updated 
            else
                error('unknown message type.')
            end
        end
    end

    function collect_data(i,ps)
        % Agents get measurements from the power grid. They update their model
        % based on the actual values from the ps structure.
        % update generators on the bus
        ig = (ps_agents(i).gen(:,C.ge.bus) == ps_agents(i).bus_id);
        ps_agents(i).gen(ig,C.ge.Pg) = ps.gen(ig,C.ge.Pg);
        ps_agents(i).gen(ig,C.ge.status) = ps.gen(ig,C.ge.status);
        % update load factor
        i_sh = (ps_agents(i).shunt(:,C.sh.bus) == ps_agents(i).bus_id);
        ps_agents(i).shunt(i_sh,C.sh.factor) = ps.shunt(i_sh,C.sh.factor);
        % update measured flows on the branches connected to the bus
        [br_id_f, ~] = find(ps_agents(i).branch(:,C.br.f) == ps_agents(i).bus_id); % find branches that this bus is on the from side
        ps_agents(i).branch(br_id_f,C.br.Pf) = ps.branch(br_id_f,C.br.Pf);
        [br_id_t, ~] = find(ps_agents(i).branch(:,C.br.t) == ps_agents(i).bus_id); % find branches that this bus is on the to side
        ps_agents(i).branch(br_id_t,C.br.Pt) = ps.branch(br_id_t,C.br.Pt);
        br_id = [br_id_f;br_id_t];
        ps_agents(i).branch(br_id,C.br.status) = ps.branch(br_id,C.br.status);
    end

    function reset_agent(i,ramp_limits)
        % reset all capacities, delta_Pg0 and delta_Pd0 for all agents 
        % generation and ramping
        ig = (ps_agents(i).gen(:,C.ge.bus) == ps_agents(i).bus_id);
        ps_agents(i).capacity.Pg = ps_agents(i).gen(ig,C.ge.Pg) .* ps_agents(i).gen(ig,C.ge.status);
        ps_agents(i).capacity.this_ramp = ramp_limits(ig) .* ps_agents(i).gen(ig,C.ge.status);
        ps_agents(i).gen(:,C.ge.ramp_rate_down) = ramp_limits;
        % loads
        ish = (ps_agents(i).shunt(:,C.sh.bus) == ps_agents(i).bus_id);
        ps_agents(i).capacity.sf = ps_agents(i).shunt(ish,C.sh.factor);
        % delta_Pg0_all and delta_Pd0_all (a model for the whole ps)
        ng = size(ps_agents(i).gen,1);
        nd = size(ps_agents(i).shunt,1);
        ps_agents(i).delta_Pg0_all = zeros(ng,1);
        ps_agents(i).delta_Pd0_all = zeros(nd,1);
        % delta_Pg and delta_sf (only for the agent -> to be implemented)
        ps_agents(i).delta_Pg = zeros(sum(ig),1);
        ps_agents(i).delta_sf = zeros(sum(ish),1);
    end

    function update_agent_data(local_nei,delta_Pd,delta_Pg)
        for j = local_nei
            bus_id = ps_agents(j).bus_id;
            % find the generation indices and update gen capacity and delta gen 
            ig = (ps_agents(j).gen(:,C.ge.bus) == bus_id);
            ps_agents(j).capacity.Pg = ps_agents(j).capacity.Pg + delta_Pg(ig);
            ps_agents(j).delta_Pg = ps_agents(j).delta_Pg + delta_Pg(ig);
            % update ramping capacity
            ps_agents(j).capacity.this_ramp = ps_agents(j).capacity.this_ramp + delta_Pg(ig);

            % find the load indices and update shunt capacity and delta load
            ish = (ps_agents(j).shunt(:,C.sh.bus) == bus_id);
            delta_sf = delta_Pd./ps_agents(j).shunt(:,C.sh.P);
            r_zero_loads = (ps_agents(j).shunt(:,C.sh.P) == 0); % correct delta_sf for loads with zero power
            delta_sf(r_zero_loads) = 0;
            ps_agents(j).capacity.sf = ps_agents(j).capacity.sf + delta_sf(ish);
            ps_agents(j).delta_sf = ps_agents(j).delta_sf + delta_sf(ish);
        end
    end

end
