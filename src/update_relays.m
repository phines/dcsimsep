function [relay,br_out_new,dt,n_over] = update_relays(ps,verbose,dt_max)
% usage: [relay,br_out_new,dt,n_over] = update_relays(ps,verbose,dt_max)

% default inputs
if nargin<2, verbose=1; end
if nargin<3, dt_max=Inf; end

% initialize outputs
C = psconstants;
relay = ps.relay;
br_out_new = [];
%dt = 0;
SMALL_EPS = 1e-12;
BIG_EPS = 1e-6;

% initialize
nr = size(ps.relay,1);
trip_t = zeros(nr,1);
excess = zeros(nr,1);
n_over = 0;

% figure out how much time until each relay trips
for r =1:nr
    switch relay(r,C.re.type)
        case C.re.oc % if this is an overcurrent relay
            % update overcurrent relays
            br_i = relay(r,C.re.branch_loc);
            bus_loc = relay(r,C.re.bus_loc);
            Imax = relay(r,C.re.setting1);
            if bus_loc==0
                Imag = max(ps.branch(br_i,C.br.Imag_f),ps.branch(br_i,C.br.Imag_t));
            elseif bus_loc==ps.branch(br_i,1)
                Imag = ps.branch(br_i,C.br.Imag_f);
            else
                Imag = ps.branch(br_i,C.br.Imag_t);
            end
            excess(r) = (Imag - Imax);
            if Imag > Imax + BIG_EPS
                n_over = n_over+1;
            end
            dist_to_threshold = relay(r,C.re.threshold) - relay(r,C.re.state_a);
            if excess(r)<=0
                trip_t(r) = Inf;
                relay(r,C.re.state_a) = relay(r,C.re.state_a) - BIG_EPS; % tick down the relay state slightly
            else
                trip_t(r) = dist_to_threshold/excess(r);
            end
        case C.re.dist
            % update distance relays
            error('not working yet');
        case C.re.uv
            % update undervoltage relays
            error('not working yet');
        otherwise
            error('Unknown relay type');
    end
end

% figure out the time step size
[dt_trip,~] = min(trip_t(trip_t>0));
%{
if dt_trip > 1e6
    dt = Inf; % This stops simulation once it gets out of update_relays 
else
    dt = min(dt_max,dt_trip);
end
%}
dt = min(dt_max,dt_trip); % to see control actions, we should only use this even if dt_trip > 1e6

% update all of the relays based on their excess
relay(:,C.re.state_a) = max(relay(:,C.re.state_a) + excess*dt + SMALL_EPS,0);

% check to see if any relays trip
if (dt_trip <= dt_max) && (dt_trip ~= Inf)
    relay_trips = find(relay(:,C.re.state_a)>=relay(:,C.re.threshold));
    % error check:
%     diff_in_trip_t = trip_t(relay_trips)-trip_t(relay_trips(1));
%     if any(abs(diff_in_trip_t)>BIG_EPS)
%         error('Something strange happened.');
%     end
    % deal with relays:
    for r = relay_trips'
        switch relay(r,C.re.type)
            case C.re.oc % overcurrent
                % trip the branch
                br_i = relay(r,C.re.branch_loc);
                br_out_new = cat(1,br_out_new,br_i);
                relay(r,C.re.state_a) = 0;

                % print something
                if verbose
                    Imag = max(ps.branch(br_i,C.br.Imag_f),ps.branch(br_i,C.br.Imag_t));
                    Imax = relay(r,C.re.setting1);
                    if Imax<BIG_EPS
                        error('Relay upper limit is set to zero');
                    end
                    fprintf(' Branch %d will trip on overcurrent in %g sec. (%g>%g)\n',br_i,dt,Imag,Imax);
                end
            otherwise
                error('Unknown relay type');
        end
    end
end
