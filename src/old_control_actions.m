function [ps, MW_lost, imbalance] = old_control_actions(ps,sub_grids,ramp_rate,it_no,opt)

% Constants
C = psconstants;
EPS = 1e-4;
% Collect some data from the system
n = size(ps.bus,1);
m = size(ps.branch,1);
F = ps.bus_i(ps.branch(:,C.br.from));
T = ps.bus_i(ps.branch(:,C.br.to));
flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,C.br.rateB);
ge_status = ps.gen(:,C.ge.status);
Pg_max = ps.gen(:,C.ge.Pmax).*ge_status + EPS;
Pg_min = ps.gen(:,C.ge.Pmin).*ge_status - EPS;
G = ps.bus_i(ps.gen(:,1));
D = ps.bus_i(ps.shunt(:,1));
Pd0 = ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.factor);
imbalance = 0;

% If we are to use the comm model do:
if opt.sim.use_comm_model
    pid = feature('getpid');
    if ispc
        comm_status_file = sprintf('C:\Temp\comm_status_%d.csv',id);
        grid_status_file = sprintf('C:\Temp\grid_status_%d.csv',pid);
    else
        comm_status_file = sprintf('/tmp/comm_status_%d.csv',pid);
        grid_status_file = sprintf('/tmp/grid_status_%d.csv',pid);
    end

    % Write out the status of each node in the network
    station_service_status = find_buses_with_power(ps,opt);
    ps.bus(:,C.bu.status) = station_service_status;
    %  output file name is grid_status_{pid}.csv
    % The comm system is only affected if there is a bi-directional link
    grid_to_comm_status = station_service_status | (~ps.bus(:,C.bu.grid_comm));
    % Write a header to the file
    fileID = fopen(grid_status_file,'w');
    fprintf(fileID,'bus,status\n');
    fclose(fileID);
    % Write the data
    dlmwrite(grid_status_file, [ps.bus(:,1) grid_to_comm_status],'-append','delimiter',',');
    % Call python comms code
    if ~all(grid_to_comm_status)
        % Figure out where the python code is
        python_location = opt.comm.python_location;
        comm_model = opt.comm.comm_model;
        % if it exists, call the comm model code
        if exist(comm_model,'file')
            systemcall = sprintf('%s %s %d %d %d',python_location,comm_model,pid,it_no,-1);
            system(systemcall);
        end
    end
    % check the status of the comm network
    % read the file /tmp/comm_status_{pid}.csv
    % bus,status
    % 1,1
    % 2,0 etc...
    if exist(comm_status_file,'file')
        data = csvread(comm_status_file,1);
        comm_status = (data(:,2)==1) | ps.bus(:,C.bu.grid_comm);
        if length(comm_status)~=n
            error('Wrong number of items in comm status file');
        end
        delete(comm_status_file);
    else
        comm_status = true(n,1);
    end
    % Read power system data (flow) from each operating node
    is_br_readable = comm_status(F) | comm_status(T);
    measured_flow = nan(m,1);
    measured_flow(is_br_readable) = flow(is_br_readable);
    measured_branch_st = true(m,1);
    measured_branch_st(is_br_readable) = ps.branch(is_br_readable,C.br.status);
else
    measured_flow = flow;
    measured_branch_st = ps.branch(:,C.br.status);
    comm_status = true(n,1);
end
comm_failed_set = find(~comm_status);

% if there are overloads in the system, try to mitigate them
if any(abs(measured_flow)>flow_max)
    % Check the mismatch
    mis_old = total_P_mismatch(ps);
    if abs(mis_old)>EPS
        error('System not balanced on entry to take_control_actions');
    end
    % Figure out the ramp rate
    ramp_dt = opt.sim.fast_ramp_mins * 60; 
    max_ramp = ramp_rate*ramp_dt;
    % Find the optimal load/gen shedding
    if opt.sim.use_mpc
        [delta_Pd, delta_Pg] = mpc_smp_solve_lp(ps,[],max_ramp,opt);
    else
        [delta_Pd, delta_Pg] = emergency_control(ps,measured_flow,measured_branch_st,max_ramp,comm_status,opt);
    end
    % If emergency control says that we should do something:
    if any(abs(delta_Pd)>EPS)
        % check to see which loads/gens can be controlled
        is_D_failed = ismember(D,comm_failed_set);
        delta_Pd(is_D_failed) = 0;
        is_G_failed = ismember(G,comm_failed_set);
        delta_Pg(is_G_failed) = 0;
        % Compute the new amount of generation
        Pg_new = ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status) + delta_Pg;
        % Error checking
%        if any( Pg_new<Pg_min | Pg_new>Pg_max )
        if any( round(Pg_new)<round(Pg_min-EPS) | round(Pg_new)>round(Pg_max+EPS) )
            error('Pg_new is out of bounds');
        end
        % Compute the new load factor
        delta_lf = delta_Pd./ps.shunt(:,C.sh.P);
        delta_lf(isnan(delta_lf)) = 0;
        lf_new = ps.shunt(:,C.sh.factor) + delta_lf;
        % Implement the results
        ps.gen(:,C.ge.P) = max(Pg_min,min(Pg_new,Pg_max)); % implement Pg
        ps.shunt(:,C.sh.factor) = max(0,min(lf_new,1));
        % Get the new mismatch
        mis_new = total_P_mismatch(ps,sub_grids);
        imbalance = imbalance + mis_new;
        % If there was an error in the balance, run rebalance again
        if abs(mis_new)>EPS
            ps = rebalance(ps,sub_grids,max_ramp,opt);
            ge_status = ps.gen(:,C.ge.status);
            Pg_max = ps.gen(:,C.ge.Pmax).*ge_status + EPS;
            Pg_min = ps.gen(:,C.ge.Pmin).*ge_status - EPS;
        end
        % Run dcpf to get subgrids
        ps = dcpf(ps,sub_grids);
        % Check to make sure things are within bounds
        Pg = ps.gen(:,C.ge.Pg);
%         if any( Pg>Pg_max | Pg<Pg_min )
        if any( round(Pg)<round(Pg_min-EPS) | round(Pg)>round(Pg_max+EPS) )
            error('Pg is out of bounds'); 
        end
    end
    Pd = ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.factor);
    Pg = ps.gen(:,C.ge.P) .* ps.gen(:,C.ge.status);
    if abs(sum(Pd)-sum(Pg)) > EPS
        error('Something is wrong.')
    else
        MW_lost = sum(Pd0) - sum(Pd);
    end
else
    MW_lost = 0;
end




