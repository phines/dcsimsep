function [out1,ge_status,d_factor] = rebalance(ps,sub_grids,ramp_limits,opt,subgrid_list,pf_model)
% very simple rebalance function
% usage: [Pg,ge_status,d_factor] = rebalance(ps,sub_grids,ramp_limits,opt,subgrid_list,pf_model)
%  or more simply:
% ps = rebalance(ps,sub_grids,ramp_limits,opt,subgrid_list,pf_model)
% 
% Inputs:
%  ps - the case data
%  sub_grids - an nx1 vector of integers indication which island each bus
%   is in
%  ramp_limits - an ngx1 vector of up/down ramp limits for the system
%  opt - global options
%  subgrid_list - the list of sub-grids that need a rebalance
%  pf_model - can be 'ac' or 'dc' (default) depending on the model used

C = psconstants;
EPS = 1e-6;
if ~isfield(ps,'bus_i');
    ps = updateps(ps);
end

% check the inputs
if nargin<2 || isempty(sub_grids)
    nodes = ps.bus(:,C.bu.id);
    br_st = (ps.branch(:,C.br.status) == 1);
    links = [ps.branch(br_st,C.br.from), ps.branch(br_st,C.br.to)];
    sub_grids = find_subgraphs(nodes,links);
end
ge_status = ps.gen(:,C.ge.status)==1;
if nargin<3 || isempty(ramp_limits), ramp_limits=ps.gen(:,C.ge.Pmax).*ge_status; end
if nargin<4, opt=psoptions; end
if nargin<5, subgrid_list = []; end
if nargin<6, pf_model = 'dc'; end % for legacy reasons

% get the loss factor
if strcmp(pf_model,'ac')
    loss_factor = 1 + opt.pf.loss_factor; % to multiply by load
elseif strcmp(pf_model,'dc')
    loss_factor = 1;
else
    error('Unknown power flow model.')
end

% collect the load data
D = ps.bus_i(ps.shunt(:,1));
Pd = ps.shunt(:,C.sh.P);
d_factor = ps.shunt(:,C.sh.factor);
% collect the generator data
G = ps.bus_i(ps.gen(:,1));
rr = ramp_limits;
Pg0 = ps.gen(:,C.ge.P).*ge_status;
Pg  = Pg0; % make a copy to modify
Pg_max = min(ps.gen(:,C.ge.Pmax),Pg+ramp_limits).*ge_status;
Pg_min = max(ps.gen(:,C.ge.Pmin),Pg-ramp_limits).*ge_status;

% check generation if needed
if opt.pf.check_Pg
    if any(Pg<Pg_min-EPS) || any(Pg>Pg_max+EPS)
        error('Generation outside of Pmin/Pmax before rebalance.');
    end
end

% find sub-grids if they are not passed in
n_sub = max(sub_grids);
if isempty(subgrid_list)
    subgrid_list = 1:n_sub; % loop over all subgrids
else
    subgrid_list = subgrid_list(:)'; % make sure it is a row vector
end

for g = subgrid_list
    bus_set = find(g==sub_grids);
    Gsub = ismember(G,bus_set) & ge_status;
    Dsub = ismember(D,bus_set) & d_factor>0;
    Pg_sub = sum(Pg(Gsub));
    Pd_sub = sum(Pd(Dsub).*d_factor(Dsub));
    % if there is no imbalance, nothing to do
    if abs(Pg_sub-Pd_sub*loss_factor)<EPS
        continue;
    end
    if opt.verbose
        if Pg_sub>Pd_sub*loss_factor
            fprintf(' Rebalance: Attempting to correct %.4f MW gen surplus in subgrid %d of %d.\n',Pg_sub-Pd_sub*loss_factor,g,n_sub);
        else
            fprintf(' Rebalance: Attempting to correct %.4f MW gen deficit in subgrid %d of %d.\n',Pd_sub*loss_factor-Pg_sub,g,n_sub);
        end
    end
    % if there are no generators in this island
    if ~any(Gsub)
        d_factor(Dsub) = 0;
        if opt.verbose
            shed_load = Pd_sub;
            fprintf(' Rebalance: No generators in subgrid %d of %d. Shed %4.2f MW of load.\n',g,n_sub,shed_load);
        end
        continue;
    end
    % if there are no loads in this island
    if ~any(Dsub) || Pd_sub<=0
        ge_status(Gsub) = 0;
        if opt.verbose
            fprintf(' Rebalance: No loads in subgrid %d of %d. Tripping generator on bus(es):\n',g,n_sub)
            all_bus = full(G(Gsub));
            all_Pg_sub = full(Pg(Gsub));
            for i = 1:length(all_bus)
                this_busID = find(ps.bus_i == all_bus(i));
                fprintf(' #%d, Pg = %.2f MW\n',this_busID,all_Pg_sub(i));
            end
        end
        Pg(Gsub) = 0;
        if ~isempty(Dsub)
            d_factor(Dsub)=0;
        end
        % shut the buses down
        ps.bus(bus_set,C.bu.status) = 0;
        continue;
    end 
    % if we want to do this the simple way. Make sure Pg_min is set to
    % zero
    if opt.sim.simple_rebalance
        if Pg_sub>Pd_sub*loss_factor
            gen_factor = Pd_sub*loss_factor/Pg_sub;
            Pg(Gsub) = Pg(Gsub) * gen_factor;
        else
            load_factor = Pg_sub/(Pd_sub*loss_factor);
            d_factor(Dsub) = d_factor(Dsub) * load_factor;
        end
        % debug
        Pg_sub = sum(Pg(Gsub));
        Pd_sub = sum(Pd(Dsub).*d_factor(Dsub));
        if abs(Pg_sub - Pd_sub*loss_factor)>EPS
            error('This shouldn''t happen');
        end
        continue
    end

    % if there is too much generation, ramp down generation
    while (Pg_sub-Pd_sub*loss_factor)>+EPS
        % figure out which generators can ramp down
        ramp_set = find(Gsub & Pg>Pg_min  & ge_status);
        if isempty(ramp_set) % break if no generators can ramp
            break
        end
        % decrease the generation that is available
        %Pg(ramp_set) = Pg0(ramp_set); % get back to the starting point
        factor = (Pd_sub*loss_factor-sum(Pg(Gsub)))/sum(rr(ramp_set));
        if factor>0
            error('Something is fishy');
        end
        Pg_old = Pg; % record old Pg for printing
        Pg(ramp_set) = min( max( Pg_min(ramp_set), Pg(ramp_set)+rr(ramp_set)*factor ), Pg0(ramp_set) );
        Pg_sub = sum(Pg(Gsub));
        if opt.verbose
            all_ramp_bus = unique(ps.gen(ramp_set,C.ge.bus));
            n = length(all_ramp_bus);
            if n > 10
                fprintf(' Rebalance: Ramping (down) generators on %d buses.\n',n)
            else
                fprintf(' Rebalance: Ramping (down) generators on bus: \n');
                for i = 1:length(all_ramp_bus)
                    idx_ramp_set = ps.gen(ramp_set,C.ge.bus) == all_ramp_bus(i);
                    idx_Pg = ramp_set(idx_ramp_set);
                    this_ramp = sum(Pg_old(idx_Pg) - Pg(idx_Pg));
                    fprintf('   #%d for %.2f MW\n',all_ramp_bus(i),this_ramp);
                end
            end
        end
    end
    % if we still have too much, trip generators
    while (Pg_sub-Pd_sub*loss_factor)>+EPS
        % figure out which generators are active
        genset = find(Gsub & ge_status & Pg>0);
        if isempty(genset) % break if no generators can be shut down
            break
        end
        % trip the smallest generator in ramp_set
        [~,i] = min(Pg(genset));
        gi = genset(i);
        if opt.verbose
            fprintf(' rebalance: Tripping generator on bus: %d. Pg=%.2f\n',ps.gen(gi,C.ge.bus),Pg(gi));
        end
        Pg(gi) = 0;
        ge_status(gi) = 0;
        Pg_sub = sum(Pg(Gsub));
    end
    if (Pg_sub-Pd_sub*loss_factor) > EPS
        if Pg_sub<EPS && Pd_sub<0
            % probably means that there are negative loads. Shut them down.
            ge_status(Gsub) = 0;
            Pg(Gsub) = 0;
            d_factor(Dsub)=0;
            % shut the buses down
            ps.bus(:,C.bu.status) = 0;
            continue
        end
        error('We should not be here');
    end
    % if there is too little generation:
    while (Pg_sub-Pd_sub*loss_factor)<-EPS
        % figure out which generators can ramp up
        ramp_set = find(Gsub & Pg<Pg_max & ge_status);
        if isempty(ramp_set) % break if no generators can ramp
            break
        end
        % increase the generation that is available
        %Pg(ramp_set) = Pg0(ramp_set).*ge_status(ramp_set); % get back to the starting point
        factor = (Pd_sub*loss_factor-sum(Pg(Gsub)))/sum(rr(ramp_set));
        if factor<0
            error('Something is fishy');
        end
        Pg_old = Pg; % record old Pg for printing
        Pg(ramp_set) = min( Pg(ramp_set) + rr(ramp_set)*factor, Pg_max(ramp_set) );
        Pg_sub = sum(Pg(Gsub));
        if opt.verbose
            all_ramp_bus = unique(ps.gen(ramp_set,C.ge.bus));
            n = length(all_ramp_bus);
            if n > 10
                fprintf(' Rebalance: Ramping (up) generators on %d buses.\n',n)
            else
                fprintf(' Rebalance: Ramping (up) generators on bus: \n');
                for i = 1:length(all_ramp_bus)
                    idx_ramp_set = ps.gen(ramp_set,C.ge.bus) == all_ramp_bus(i);
                    idx_Pg = ramp_set(idx_ramp_set);
                    this_ramp = sum(Pg(idx_Pg) - Pg_old(idx_Pg));
                    fprintf('   #%d for %.2f MW\n',all_ramp_bus(i),this_ramp);
                end
            end
        end
    end
    % if we still have too little, do load shedding
    Pd_sub0 = Pd_sub; 
    if Pd_sub*loss_factor > (Pg_sub + EPS)
        factor = Pg_sub/(Pd_sub*loss_factor);
        d_factor(Dsub) = factor * d_factor(Dsub);
        Pd_sub = sum(d_factor(Dsub).*Pd(Dsub));
        if opt.verbose
            shed_load = Pd_sub0 - Pd_sub;
            fprintf(' Rebalance: Shed %4.2f MW of load in subgrid %d of %d.\n',shed_load,g,n_sub);
        end
    end
    % check for a balance
    if abs(Pg_sub-Pd_sub*loss_factor)>EPS
        error('Could not find a balance');
    end
end
% debug check
if opt.pf.check_Pg
    Pg_min = Pg_min.*ge_status;
    Pg_max = Pg_max.*ge_status;
    if any(Pg<Pg_min-EPS) || any(Pg>Pg_max+EPS)
        error('Generation outside of Pmin/Pmax after rebalance.');
    end
end

% Check the number of outputs
if nargout==1
    % In single output version, out1=ps
    ps.shunt(:,C.sh.factor) = d_factor;
    ps.gen(:,C.ge.status) = ge_status;
    ps.gen(:,C.ge.P) = Pg;
    out1 = ps;
else
    % in multiple output version, out1=Pg
    out1 = Pg;
end

end