function event = make_event(branch_outages,bus_outages,t_outages)
% make an event matrix for a set of outages that occur at t_outages

C = psconstants;

n_br_outs = length(br_outages);
n_bu_outs = length(bus_outages);
br_rows = (1:n_br_outs) + 1;
bu_rows = (1:n_bu_outs) + 1 + n_br_outs;

% init the event matrix
event = zeros(n_br_outs+n_bu_outs+2,C.ev.cols);
% branch outages
event(br_rows,C.ev.time) = 1;
event(br_rows,C.ev.type) = C.ev.trip_branch;
event(br_rows,C.ev.branch_loc) = br_outages;
% bus outages
event(bu_rows,C.ev.time) = 1;
event(bu_rows,C.ev.type) = C.ev.trip_bus;
event(bu_rows,C.ev.bus_loc) = bus_outages;

% finish
event(end,C.ev.type) = C.ev.finish;
event(end,C.ev.time) = 3600; % one hour
