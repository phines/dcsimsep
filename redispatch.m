function [out1,ge_status,d_factor] = redispatch(ps,sub_grids,ramp_limits,verbose,opt)
% very simple redispatch function
% usage: [Pg,ge_status,d_factor] = redispatch(ps,sub_grids,ramp_limits,verbose)
%  or more simply:
% ps = redispatch(ps,sub_grids,ramp_limits,verbose)
% 
% Inputs:
%  ps - the case data
%  sub_grids - an nx1 vector of integers indication which island each bus
%   is in
%  ramp_limits - an mx1 vector of up/down ramp limits for the system
%  verbose - flag for output

C = psconstants;
EPS = 1e-6;
n = size(ps.bus,1);
if ~isfield(ps,'bus_i');
    ps = updateps(ps);
end

% check the inputs
if nargin<2||isempty(sub_grids), sub_grids=ones(n,1); end
ge_status = ps.gen(:,C.ge.status)==1;
if nargin<3||isempty(ramp_limits), ramp_limits=ps.gen(:,C.ge.Pmax).*ge_status; end
if nargin<4, verbose=true; end
if nargin<5, opt=psoptions; end

% collect the load data
D = ps.bus_i(ps.shunt(:,1));
Pd = ps.shunt(:,C.sh.P);
d_factor = ps.shunt(:,C.sh.factor);
% collect the generator data
G = ps.bus_i(ps.gen(:,1));
rr = ramp_limits;
%ge_status = ps.gen(:,C.ge.status)==1;
Pg0 = ps.gen(:,C.ge.P).*ge_status;
Pg  = Pg0;
Pg_max = min(ps.gen(:,C.ge.Pmax),Pg+ramp_limits).*ge_status;
Pg_min = max(ps.gen(:,C.ge.Pmin),Pg-ramp_limits).*ge_status;
% if any(Pg<Pg_min-EPS) || any(Pg>Pg_max+EPS)
if any(round(Pg)<round(Pg_min-EPS)) || any(round(Pg)>round(Pg_max+EPS))
    keyboard
    error('Generation outside of Pmin/Pmax');
end

grid_list = unique(sub_grids)';
n_sub = length(grid_list);
for g = grid_list
    bus_set = find(g==sub_grids);
    Gsub = ismember(G,bus_set) & ge_status;
    Dsub = ismember(D,bus_set) & d_factor>0;
    Pg_sub = sum(Pg(Gsub));
    Pd_sub = sum(Pd(Dsub).*d_factor(Dsub));
    % if there is no imbalance, nothing to do
    if abs(Pg_sub-Pd_sub)<EPS
        continue;
    end
    if verbose
        if Pg_sub>Pd_sub
            fprintf('Attempting to correct %.4f MW gen surplus in subgrid %d of %d\n',Pg_sub-Pd_sub,g,n_sub);
        else
            fprintf('Attempting to correct %.4f MW gen deficit in subgrid %d of %d\n',Pd_sub-Pg_sub,g,n_sub);
        end
    end
    % if there are no generators in this island
    if ~any(Gsub)
        d_factor(Dsub) = 0;
        continue;
    end
    % if there are no loads in this island
    if ~any(Dsub) || Pd_sub<0
        ge_status(Gsub) = 0;
        Pg(Gsub) = 0;
        if ~isempty(Dsub)
            d_factor(Dsub)=0;
        end
        % shut the buses down
        ps.bus(bus_set,C.bu.status) = 0;
        continue;
    end
    % if we want to do this the simple way
    if opt.sim.simple_redispatch
        if Pg_sub>Pd_sub
            gen_factor = Pd_sub/Pg_sub;
            Pg(Gsub) = Pg(Gsub) * gen_factor;
        else
            load_factor = Pg_sub/Pd_sub;
            d_factor(Dsub) = d_factor(Dsub) * load_factor;
        end
        % debug
        Pg_sub = sum(Pg(Gsub));
        Pd_sub = sum(Pd(Dsub).*d_factor(Dsub));
        if abs(Pg_sub - Pd_sub)>EPS
            error('This shouldn''t happen');
        end
        continue
    end

    % if there is too much generation, ramp down generation
    while (Pg_sub-Pd_sub)>+EPS
        % figure out which generators can ramp down
        ramp_set = find(Gsub & Pg>Pg_min  & ge_status);
        if isempty(ramp_set) % break if no generators can ramp
            break
        end
        % decrease the generation that is available
        %Pg(ramp_set) = Pg0(ramp_set); % get back to the starting point
        factor = (Pd_sub-sum(Pg(Gsub)))/sum(rr(ramp_set));
        if factor>0
            error('Something is fishy');
        end
        if verbose
            fprintf(' Ramping (down) generators on bus(es): ');
            fprintf('#%d ',unique(ps.gen(ramp_set,C.ge.bus)));
            fprintf('\n');
        end
        Pg(ramp_set) = min( max( Pg_min(ramp_set), Pg(ramp_set)+rr(ramp_set)*factor ), Pg0(ramp_set) );
        Pg_sub = sum(Pg(Gsub));
    end
    % if we still have too much, trip generators
    while (Pg_sub-Pd_sub)>+EPS
        % figure out which generators are active
        genset = find(Gsub & ge_status & Pg>0);
        if isempty(genset) % break if no generators can be shut down
            break
        end
        % trip the smallest generator in ramp_set
        [~,i] = min(Pg(genset));
        gi = genset(i);
        if verbose
            fprintf(' Tripping generator on bus: %03d. Pg=%.2f\n',ps.gen(gi,C.ge.bus),Pg(gi));
        end
        Pg(gi) = 0;
        ge_status(gi) = 0;
        Pg_sub = sum(Pg(Gsub));
    end
    if (Pg_sub-Pd_sub) > EPS
        if Pg_sub<EPS && Pd_sub<0
            % probably means that there are negative loads. Shut them down.
            ge_status(Gsub) = 0;
            Pg(Gsub) = 0;
            d_factor(Dsub)=0;
            % shut the buses down
            ps.bus(:,C.bu.status) = 0;
            continue;
        end
        error('We should not be here');
    end
    % if there is too little generation:
    while (Pg_sub-Pd_sub)<-EPS
        % figure out which generators can ramp up
        ramp_set = find(Gsub & Pg<Pg_max & ge_status);
        if isempty(ramp_set) % break if no generators can ramp
            break
        end
        % increase the generation that is available
        %Pg(ramp_set) = Pg0(ramp_set).*ge_status(ramp_set); % get back to the starting point
        factor = (Pd_sub-sum(Pg(Gsub)))/sum(rr(ramp_set));
        if factor<0
            error('Something is fishy');
        end
        if verbose
            fprintf(' Ramping (up) generators on bus: ');
            fprintf('#%d ',unique(ps.gen(ramp_set,C.ge.bus)));
            fprintf('\n');
        end
        Pg(ramp_set) = min( Pg(ramp_set) + rr(ramp_set)*factor, Pg_max(ramp_set) );
        Pg_sub = sum(Pg(Gsub));
    end
    % if we still have too little, do load shedding
    Pd_sub0 = Pd_sub;
    if Pd_sub > (Pg_sub + EPS)
        factor = Pg_sub/Pd_sub;
        d_factor(Dsub) = factor * d_factor(Dsub);
        Pd_sub = sum(d_factor(Dsub).*Pd(Dsub));
        if verbose
            shed_load = Pd_sub0 - Pd_sub;
            fprintf(' Shed %6.2f MW of load\n',shed_load);
        end
    end
    % check for a balance
    if abs(Pg_sub-Pd_sub)>EPS
        error('Could not find a balance');
    end
end

% debug check
Pg_max = min(ps.gen(:,C.ge.Pmax),Pg+ramp_limits).*ge_status;
Pg_min = max(ps.gen(:,C.ge.Pmin),Pg-ramp_limits).*ge_status;
% if any(Pg<Pg_min-EPS) || any(Pg>Pg_max+EPS)
if any(round(Pg)<round(Pg_min-EPS)) || any(round(Pg)>round(Pg_max+EPS))
    error('Generation outside of Pmin/Pmax');
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

return;
    
    
