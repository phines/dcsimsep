function relay = relay_settings(ps,distance,overcurrent,undervoltage)
% produce rough relay settings for ps

if nargin<2, distance = true;     end
if nargin<3, overcurrent = true;  end
if nargin<4, undervoltage = true; end

% extract data
C = psconstants;
bus_nos = ps.bus(:,1);
%n = size(ps.bus,1);
m = size(ps.branch,1);
%ng = size(ps.gen);
%ns = size(ps.shunt);
%nbr = size(ps.branch);
F = ps.bus_i(ps.branch(:,1));
%T = ps.bus_i(ps.branch(:,2));

%% prepare output
relay = [];

%% add simple zone 3 distance relays
if distance
    %relay_dist = zeros(nbr*2,C.re.cols);
end

%% add simple under-voltage relays
if undervoltage
    %relay_uv = zeros(n,C.re.cols);
end

%% add simple time overcurrent relays
if overcurrent
    relay_oc_f = zeros(m,C.re.cols);
    %relay_oc_t = zeros(m,C.re.cols);
    
    % typ
    relay_oc_f(:,C.re.type) = C.re.oc;
    %relay_oc_t(:,C.re.type) = C.re.oc;
    
    % relay locations
    relay_oc_f(:,C.re.branch_loc) = (1:m)';
    relay_oc_f(:,C.re.bus_loc)    = bus_nos(F);
    %relay_oc_t(:,C.re.branch_loc) = (1:m)';
    %relay_oc_t(:,C.re.bus_loc)    = bus_nos(T);
    
    % settings
    Imax = ps.branch(:,C.br.rateB) / ps.baseMVA;
    overload_max = Imax * 1.5 * 5; % 5 seconds at a 50% overload produces a trip
    relay_oc_f(:,C.re.setting1)  = Imax;
    relay_oc_f(:,C.re.threshold) = overload_max;
    %relay_oc_t(:,C.re.setting1)  = Imax;
    %relay_oc_t(:,C.re.threshold) = overload_max;
    
    %relay = [relay;
    %    relay_oc_f;
    %    relay_oc_t];
    relay = [relay;relay_oc_f];
end

% done
return
