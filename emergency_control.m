function [delta_Pd,delta_Pg] = emergency_control(ps,sub_grids,ramp_limits,measured_flow,comm_status,verbose)
% An emergency control function to adjust load, based on current measured
%  flows
% usage: [delta_df,delta_Pg] = emergency_opf(ps,sub_grids,ramp_limits,measured_flow,verbose)
% inputs:
%  ps - power systems data
%  sub_grids - indicates which island/subgrid each bus is in.
%  ramp_limits - the maximum amount (in MW) that each generator can ramp in
%   this time period
%  measured_flow - the amount of flow (in MW/MVA) on each transmission line
%  comm_status - a binary vectory indicating the status of each bus'es comm
%   system
%  verbose - a binary indicating whether to print stuff about what is going
%   on.
%  
% outputs:
%  delta_Pd - changes to amount of load at each "shunt"
%  delta_Pg - changes to the generator output at each generator

C = psconstants;
EPS = 1e-6;
n = size(ps.bus,1);

% check the inputs
if nargin<4, error('need at least 4 inputs'); end
if nargin<5, comm_status = true(n,1); end
if nargin<6, verbose = true; end

% size of the system
n = size(ps.bus);
m = size(ps.branch);

% collect load data
D = ps.bus_i(ps.shunt(:,1));
nd = length(D);
Pd0 = ps.shunt(:,C.sh.P);
%d_factor = ps.shunt(:,C.sh.factor);
%d_df_min = max(-d_factor,1);
%d_df_max = min(1 - d_factor,0);
% collect generator data
G = ps.bus_i(ps.gen(:,1));
ng = length(G);
rr = ramp_limits; % MW per second I think...
ge_status = ps.gen(:,C.ge.status)==1;
Pg0 = ps.gen(:,C.ge.P).*ge_status;
Pg  = Pg0;
Pmax = ps.gen(:,C.ge.Pmax);
Pmin = zeros(ng,1);
dPg_max = max(min(Pmax-Pg,rr),0);
dPg_min = max(-Pg,-rr);
if any(Pg<Pmin-EPS) || any(Pg>Pmax+EPS)
    error('Generation outside of Pmin/Pmax');
end
% initialize the outputs
delta_Pg = zeros(ng,1);
delta_Pd = zeros(nd,1);
return

% find the overloads
flow_max = ps.branch(:,C.br.rateB)/ps.baseMVA;
overloads = find(measured_flow > flow_max);
ex_post_flow = measured_flow;

% build the ptdf matrix
ptdf = get_ptdf(ps);

% look at each overload and find the best way to reduce it
% continue to shed load and generation until the overloads are gone
for o = overloads'
    ptdf_o = ptdf(o,:);
    % how much overload do we need to eliviate?
    overload_amount = ex_post_flow(o) - flow_max(o);
    
    % find the best generator to reduce
    ptdf_g = ptdf_o(G);
    [dFlow_dPg,gi] = max(ptdf_g);
    % how much dPg is needed to remove the overload
    dPg_needed = - dFlow_dPg * overload_amount
        
    % find the best load to shed
    ptdf_d = ptdf_o(D);
    [dFlow_dPd,di] = min(ptdf_d);
    
    % determine how much to shed
    
    
    % shed them
end

return;
    
    
