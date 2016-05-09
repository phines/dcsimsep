function msg_ext = prepare_msg_ext(ps_agent,from_bus_id)
% prepare the extended message for the sending bus. The message includes
% the list of tripped lines
C = psconstants;
if nargin < 2
    msg_ext.bus_id = ps_agent.bus_id; % the bus id of the sending bus
else
    msg_ext.bus_id = from_bus_id; % this is used when information for 
                                  % another bus is gathered, e.g., in
                                  % the function share_information
end
[br_id, ~] = find(ps_agent.branch(:,[C.br.f,C.br.t]) == msg_ext.bus_id); % find branches connected to this bus
% save status of all the above branches in the msg
msg_ext.branch.status = [br_id, ps_agent.branch(br_id,C.br.status)];