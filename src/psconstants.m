function C_out = psconstants
% constants for a power system
%  These help us figure out where the data is in the data structure

% making C persistent allows us to call this function very quickly:
persistent C;

if isempty(C)
	%% constants used to indicate the type of generator control used
	C.PQ  = 1; % when a generator output is constant P/Q
	C.PV  = 2; % when a generator is set to a voltage set point
	C.REF = 3; % use for an ideal voltage source
    C.ISO = 4; % isolated bus
    C.PF  = 5; % participation factor gen/bus
	C.INF = Inf; % use to indicate an infinite bus
    
	%% constants used for devices that may be controllable:
	C.NO_CONTROL = 1;
	C.BINARY     = 2;
	C.CONTINUOUS = 3;
	
	%% constants for branches and circuit breaker status
	C.OPEN   = 0;
	C.CLOSED = 1;
	C.F_OPEN = 2; % from end open, to end closed
	C.T_OPEN = 3; % to end open, from end closed
	
	%% other constants
	C.EMTPY = -99999999; % used to indicate empty or un-filled values
    C.DEFAULT_VALUE = 1000; % default value of a load ($/MVA)
	
	%% bus data (same as MATPOWER, though loads are not used---remain zeros)
	C.bu.id     = 1;
    C.bu.type   = 2;  % use with caution---can be confusing
    C.bu.Pd     = 3;  % !!! DO NOT USE, calculated from shunt !!!
    C.bu.Qd     = 4;  % !!! DO NOT USE, calculated from shunt !!!
    C.bu.Gs     = 5;  % !!! DO NOT USE, calculated from shunt !!!
    C.bu.Bs     = 6;  % !!! DO NOT USE, calculated from shunt !!!
    C.bu.area   = 7;  % The area number for this bus
	C.bu.Vmag   = 8;  % rms magnitude, p.u.
	C.bu.Vang   = 9;  % bus angle in degrees
	C.bu.baseKV = 10; % base kV (voltage at 1.0 pu)
    C.bu.zone   = 11; % zone
	C.bu.Vmax   = 12; % minimum recommended voltage for this bus 
	C.bu.Vmin   = 13; % maximum recommended voltage for this bus
    C.bu.lam_P  = 14; % sensitivity to injection
    C.bu.lam_Q  = 15; % sensitivity to injection
    C.bu.mu_Vmax = 16; % sensitivity to changes in limits
    C.bu.mu_Vmin = 17; % sensitivity to changes in limits
	C.bu.locX   = 18;  % location of the bus in the x dimension
	C.bu.locY   = 19;  % location of the bus in the y dimension
    C.bu.comm_status = 20; % this is a binary value indicating if the comm system is working at this node
    C.bu.status = 21; % this is a binary value indicating if bus has station service power
    C.bu.power_from_sh = 22; % indicates which shunt this bus gets its station service power from
    C.bu.grid_comm = 23; % Binary indcator that tells us whether there exists a bi-directional grid-comm link
    % synonyms
    C.bu.locs   = [C.bu.locX C.bu.locY];
    
    % column names
    C.bu.col_names = {'ID','type','Pd','Qd','Gs','Bs','area','Vmag','Vang','basekV','zone','Vmax','Vmin',...
        'lam_P','lam_Q','mu_Vx','mu_Vn','locX','locY','comm_status','bus_status','power_from_sh'};

	C.bu.cols = 22; % minimum number of columns
	C.bus = C.bu; % allows us to use C.bus or C.bu
    % default kv
    C.bu.baseKV_default = 230;

	%% branch data (columns 1-17 are the same as MATPOWER)
	C.br.from   = 1;  % "from" end bus number
	C.br.to     = 2;  % "to" end bus number
	C.br.R      = 3;  % resistance (p.u.)
	C.br.X      = 4;  % line reactance (p.u.)
	C.br.B      = 5;  % line charging for the whole line (not B/2) (p.u.)
	C.br.rateA  = 6;  % normal line rating
	C.br.rateB  = 7;  % short-term line rating
	C.br.rateC  = 8;  % emergency line rating 
	C.br.tap    = 9;  % tap ratio for transformer
	C.br.shift  = 10; % phase shift caused by transformer
	C.br.status = 11; % branch status (1, 0, C.F_OPEN, or C.T_OPEN)
    C.br.Pf     = 12; % branch from end real power flow
    C.br.Qf     = 13; % branch from end reactive power flow
    C.br.Pt     = 14; % branch to end real power flow
    C.br.Qt     = 15; % branch to end reactive power flow
    C.br.mu_f   = 16; % the sensitivity to changes to the from-end flow limit
    C.br.mu_t   = 17; % the sensitivity to changes to the to-end flow limit
	C.br.Imag_f = 18; % rms current (p.u.) on the from end of the branch
	C.br.Imag_t = 19; % rms current (p.u.) on the to end of the branch
	C.br.switchable = 20; % used to indicate that the branch is switchable
    C.br.fail_rate = 21;  % failure rate of the branch (lambda - outages per year)
    C.br.type   = 22; % branch type;
    C.br.id     = 23; % branch id
    C.br.col_names = {'from','to','R','X','B','rateA','rateB','rateC','tap','shift','status',...
        'Pf','Qf','Pt','Qt','mu_f','mu_t','Imag_f','Imag_t','switchable','fail_rate','type','id'};

	C.br.cols = 23;   % min no. of cols for branch

	% short cuts:
	C.br.f = C.br.from;
	C.br.t = C.br.to;
	C.branch  = C.br;
    C.br.Imag = (C.br.Imag_f:C.br.Imag_t);
    C.br.rates = (C.br.rateA:C.br.rateC);
    % branch type constants
    C.br.type_tline = 0;
    C.br.type_twownd = 1; % two winding transformer
    C.br.type_threewnd = 2; % one branch of a three-winding transformer
    C.br.type_redtwownd = 3; % ward reduced from two winding transformer
    C.br.type_redthreewnd = 4; % ward reduced from three winding transformer

	%% generator data (cols 1-14 same as MATPOWER)
	C.ge.bus    = 1; % the location of the generator (bus number)
	C.ge.Pg     = 2; % the actual output of the generator
	C.ge.Qg     = 3; % the actual reactive output of the generator
	C.ge.Qmax   = 4; % max reactive power output
	C.ge.Qmin   = 5; % min reactive power output
	C.ge.Vsp    = 6; % set point for the generator if PV or REF
    C.ge.mBase  = 7;% machine MVA base
	C.ge.status = 8; % the generator status
   	C.ge.Pmax   = 9; % max real power output
    C.ge.Pmin   = 10; % min real power output
    C.ge.mu_Pmax = 11; % sensitivity
    C.ge.mu_Pmin = 12; % sensitivity
    C.ge.mu_Qmax = 13; % sensitivity
    C.ge.mu_Qmin = 14; % sensitivity
	C.ge.type   = 15;  % PV, PQ, PF, or REF
	C.ge.cost   = 16;  % the marginal cost ($/MW) of power from this generator
    C.ge.part_fact = 17; % the participation factor of the generator (used only if type is PF or REF);
    C.ge.ramp_rate_up = 18; % ramp rate up
    C.ge.ramp_rate_down = 19; % down ramp rate
    C.ge.col_names = {'bus','Pg','Qg','Qmax','Qmin','Vsp','mBase','status','Pmax','Pmin',...
        'mu_Px','mu_Pn','mu_Qx','mu_Qn','type','cost','part_fact','RRU','RRD','id'};

	C.ge.cols = 20; % min no. of cols
    % synonyms
    C.gen  = C.ge;
    C.ge.P = C.ge.Pg;
    C.ge.Q = C.ge.Qg;
    
    %% values for electromechanical machine dynamics
    % The notation and per unit methods are based on Bergen & Vittal.
    C.ma.gen  = 1;  % generator number (should be one entry per generator, sequential)
    C.ma.r    = 2;  % coil-winding AC resistance
    C.ma.Xd   = 3;  % synchronous reactance (p.u.)
    C.ma.Xdp  = 4;  % transient reactance (p.u.)
    C.ma.Xdpp = 5;  % sub-transient reactance (p.u.)
    C.ma.Xq   = 6;  % Q-axis synchronous reactance (pu)
    C.ma.Xqp  = 7;  % Q-axis transient reactance (pu)
    C.ma.Xqpp = 8;  % Q-axis sub-transient reactance (pu)
    C.ma.D    = 9;  % damping constant (units?)
    C.ma.M    = 10; % machine starting time (essentially inertia. pu using method in Bergen)
    % state variables
    C.ma.Ea   = 11; % the open circuit voltage (magnitude) of the machine (field-induced voltage)
    C.ma.Eap  = 12; % the transient internal voltage of the machine. (not always needed).
    C.ma.Pm   = 13; % the mechancial power of the machine, after accounting for steady state damping (friction)
    C.ma.Pm0  = 13; % also known as P_m^0 in Bergen. (system base)
    C.ma.delta = 14; % rotor angle, with respect to the system reference angle
    C.ma.delta_dot = 15; % defined such that steady state rotor speed is zero. omega = d(delta)/dt
    
    % NOTE THAT WE HAVE NOT YET DEFINED VARIABLES FOR THE GOVERNOR OR
    % EXCITER
    
    C.ma.cols = 15;
    % synonyms
    C.ma.Ea_mag = C.ma.Ea;
    C.mac = C.ma;

	%% shunt elements (capacitors, loads, etc.)
	% A shunt essentially has a constant complex power and constant complex impedance elements
	%  in parallel between a bus and ground. "factor" or "status" can be used to scale the elements
	C.sh.bus    = 1; % the location (bus number) of the shunt element
	C.sh.P      = 2; % real power consumption (MW) at 1.0 p.u. voltage
	C.sh.Q      = 3; % reactive power consumption (MVAr) at 1.0 p.u. voltage
	C.sh.frac_S = 4; % fraction of the shunt (P/Q) that is constant complex power
	C.sh.frac_Z = 5; % fraction of the shunt (P/Q) that is constant complex impedance
                     % note that the remainder is constant current (if supported)
    C.sh.frac_Y = 5; % frac_Y is a synonym for frac_Z
	C.sh.status = 6; % indicates if the shunt is on or off
	C.sh.factor = C.sh.status; % status and factor have the same effect
	% use status with binary control or factor with continuous control
	C.sh.type   = 7; % indicates how this shunt is (or is not) controlled
	C.sh.value  = 8; % the relative value of this load ($/MVA), where MVA is abs(P+jQ)
	C.sh.frac_E = 9; % fraction of the shunt that should use an expontential load model P*Vmag^gamma;
	C.sh.gamma  = 10; % the exponent gamma in the equation above.
    C.sh.near_gen = 11; % the nearest generator to get frequency signal
    C.sh.id     = 12; % unique id number for this shunt element.
	C.sh.cols   = 12; % minimum number of columns
    % synonyms:   
	C.shunt     = C.sh; % synonym
    C.sh.col_names = {'bus','P','Q','frac_S','frac_Z','status','type','value','frac_E','gamma','near_genID','id'};

    %% event-related definitions

    % event types 
    C.event.start  = 0;
    C.event.end    = 1;
    C.event.finish = C.event.end; % same as end
    C.event.fault  = 2; % three phase to ground fault by default
    C.event.clear_fault = 3; % clear a fault
    C.event.trip_branch  = 4;
    C.event.close_branch = 5;
    C.event.trip_bus     = 6;
    
    % columns 
    C.event.time       = 1; % event time in seconds
    C.event.type       = 2; % type (start and finish are required)
    C.event.bus_loc    = 3; % the bus at which the event takes place
    C.event.branch_loc = 4; % the branch at which the event takes place
    C.event.gen_loc    = 5; % the gen at which the event takes place
    C.event.shunt_loc  = 6; % the shunt at which the event takes place
    C.event.cols = 6;
    % synonym
    C.ev = C.event;
    
    %% relay matrix defines relays
    % relay types:
    %  oc: overcurrent
    %  dist: distance
    %  uv: undervoltage
    %  ufls: underfrequency load shedding
    C.relay.oc   = 1;
    C.relay.dist = 2;
    C.relay.uv   = 3;
    C.relay.ufls = 4;
    % relay columns
    C.relay.type       = 1; % 
    C.relay.setting1   = 2; % setting 1
    C.relay.setting2   = 3; % setting 2
    C.relay.threshold  = 4; % the relay threshold. When >state, relay trips
    C.relay.state_a    = 5; % the analog state of the relay
    C.relay.state      = 5; % the analog state of the relay
    C.relay.state_b    = 6; % the binary state of the relay (open(0)/closed(1));
    C.relay.bus_loc    = 7; % the bus location for the relay
    C.relay.branch_loc = 8; % the branch location for the relay
    C.relay.gen_loc    = 9; % the gen location for the relay
    C.relay.shunt_loc  = 10; % the shunt location for the relay
    C.relay.cols = 9;
    % syn
    C.re = C.relay;
end

C_out = C;
