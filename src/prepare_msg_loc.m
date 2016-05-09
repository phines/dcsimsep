function msg_loc = prepare_msg_loc(ps_agent,from_bus_id)
% prepare the local message for the sending bus. The message includes load
% and generation on the sending bus together with flows and status of all 
% lines connected to it.
C = psconstants;
if nargin < 2
    msg_loc.bus_id = ps_agent.bus_id; % the bus id of the sending bus
else
    msg_loc.bus_id = from_bus_id; % this is used when information for 
                                  % another bus is gathered, e.g., in
                                  % the function share_information
end

% save flow for branches that this bus is on the from side 
[br_id_f, ~] = find(ps_agent.branch(:,C.br.f) == msg_loc.bus_id); % find branches that this bus is on the from side
msg_loc.branch.Pf = [br_id_f, ps_agent.branch(br_id_f,C.br.Pf)]; % save branch active flows 

% do the same for branches that this bus is on the to side 
[br_id_t, ~] = find(ps_agent.branch(:,C.br.t) == msg_loc.bus_id); % find branches that this bus is on the to side
msg_loc.branch.Pt = [br_id_t, ps_agent.branch(br_id_t,C.br.Pt)]; % save branch active flows 

% save status of all the above branches in the msg
br_id = [br_id_f;br_id_t];
msg_loc.branch.status = [br_id, ps_agent.branch(br_id,C.br.status)];

% save the load factor of the bus in the msg
i_sh = find(ps_agent.shunt(:,C.sh.bus) == msg_loc.bus_id);
msg_loc.shunt.factor = [i_sh, ps_agent.shunt(i_sh,C.sh.factor)];

% save the generation on the bus in the msg
ig = find(ps_agent.gen(:,C.ge.bus) == msg_loc.bus_id);
msg_loc.gen.Pg = [ig, ps_agent.gen(ig,C.ge.Pg)];
msg_loc.gen.status = [ig, ps_agent.gen(ig,C.ge.status)];




