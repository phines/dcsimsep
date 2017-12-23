function C_out = mpConstants
% usage: C_out = mpConstants

persistent C;

if isempty(C)
  % bus constants:
  [C.PQ, C.PV, C.REF, C.NONE, C.BUS_I, C.BUS_TYPE, C.PD, C.QD, C.GS, C.BS, C.BUS_AREA, C.VM, ...
      C.VA, C.BASE_KV, C.ZONE, C.VMAX, C.VMIN, C.LAM_P, C.LAM_Q, C.MU_VMAX, C.MU_VMIN] = idx_bus;
  C.MAC = 4; % for machine bus type
  C.LOC_X = 18;
  C.LOC_Y = 19;
  C.LMP = C.LAM_P;
  C.BUS_TIME_STAMP = 20;
  C.N_BUS_COLS = 20;
  
  % gen constants:
  [C.GEN_BUS, C.PG, C.QG, C.QMAX, C.QMIN, C.VG, C.MBASE, ...
      C.GEN_STATUS, C.PMAX, C.PMIN, C.MU_PMAX, C.MU_PMIN, C.MU_QMAX, C.MU_QMIN] = idx_gen;
  C.RUR = 15;
  C.RDR = 16;
  C.RUC = 17;
  C.RDC = 18;
  C.GEN_TIME_STAMP = 19;
  C.N_GEN_COLS = 19;
  
  % machine constants:
  C.MAC_BUS = 1;
  C.MAC_MVABASE = 2;
  C.MAC_KVBASE = 3;
  C.MAC_XD = 4;
  C.MAC_D = 5;
  C.MAC_H = 6;
  C.MAC_DPU = 7;
  C.MAC_MPU = 8;
  C.MAC_GOV = 9;
  C.MAC_RR = 10;
  C.GOV_DP_DT = 11;
  C.N_MAC_COLS = 11;
  
  % branch constants:
  [C.F_BUS, C.T_BUS, C.BR_R, C.BR_X, C.BR_B, C.RATE_A, C.RATE_B, ...
      C.RATE_C, C.TAP, C.SHIFT, C.BR_STATUS, C.PF, C.QF, C.PT, C.QT, C.MU_SF, C.MU_ST] = idx_brch;
  C.IMAG_F = 18;
  C.IMAG_T = 19;
  C.VMAG_F = 20;
  C.VMAG_T = 21;
  C.THETA_FT = 22;
  C.IMAG_PREV = 23;
  % relay related stuff:
  C.OC_PROT_LIM   = 24;
  C.OC_PROT_STATE = 25;
  C.DIST_PROT_LIM = 26;
  C.DIST_PROT_STATE = 27;
  C.BR_TIME_STAMP = 28;
  C.N_BRANCH_COLS = 28;

  % load constants:
  C.L_BUS = 1;
  C.L_VALUE = 2;
  C.L_PD = 3;
  C.L_QD = 4;
  C.L_GS = 5;
  C.L_BS = 6;
  C.L_IRE = 7;
  C.L_IIM = 8;
  C.L_STATUS = 9; % this variable can be used as a binary variable 
  C.L_FACTOR = 9; % or a continuous one to indicate the extent to which the load has
                  % been reduced from its full value
  C.L_TIME_STAMP = 10;                  

  C.N_LOAD_COLS  = 10;
                  
  % other constants
  C.EMPTY = -999999999;
  C.EPS = 1e-3;
  C.SMALL = .01;
  C.j = sqrt(-1);
  C.SIG_VMAG = .1;
  C.SIG_VANG = 10;
  
  % for mpc
  C.LOAD_QUAD_COST = 10;
  C.GEN_QUAD_COST = 10;
  C.VMIN_ST = .75;
  C.VMAX_ST = 1.25;
  C.DISC_RATE = .1;
  C.SAFETY_MARGIN = .01;
  C.INCLUDE_INORM = .8; % include branches with this current or higher
  C.MAX_T = 10;
  
  % for update
  C.DEFAULT_L_VALUE = 1000;
end

C_out = C;

function [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen
%IDX_GEN   Defines constants for named column indices to gen matrix.
%
%   [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
%   PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pg = gen(4, PG);   % get the real power output of generator 4
%    gen(:, PMIN) = 0;  % set to zero the minimum real power limit of all gens
% 
%   The index, name and meaning of each column of the gen matrix is given
%   below:
%
%   columns 1-03 must be included in input matrix (in case file)
%    1  GEN_BUS     bus number
%    2  PG          Pg, real power output (MW)
%    3  QG          Qg, reactive power output (MVAr)
%    4  QMAX        Qmax, maximum reactive power output (MVAr)
%    5  QMIN        Qmin, minimum reactive power output (MVAr)
%    6  VG          Vg, voltage magnitude setpoint (p.u.)
%    7  MBASE       mBase, total MVA base of machine, defaults to baseMVA
%    8  GEN_STATUS  status, 1 - in service, 0 - out of service
%    9  PMAX        Pmax, maximum real power output (MW)
%    10 PMIN        Pmin, minimum real power output (MW)
%   
%   columns 11-14 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    11 MU_PMAX     Kuhn-Tucker multiplier on upper Pg limit (u/MW)
%    12 MU_PMIN     Kuhn-Tucker multiplier on lower Pg limit (u/MW)
%    13 MU_QMAX     Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
%    14 MU_QMIN     Kuhn-Tucker multipcd lier on lower Qg limit (u/MVAr)

%   MATPOWER
%   $Id: idx_gen.m,v 1.6 2004/12/15 22:46:32 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
GEN_BUS     = 1;    %% bus number
PG          = 2;    %% Pg, real power output (MW)
QG          = 3;    %% Qg, reactive power output (MVAr)
QMAX        = 4;    %% Qmax, maximum reactive power output (MVAr)
QMIN        = 5;    %% Qmin, minimum reactive power output (MVAr)
VG          = 6;    %% Vg, voltage magnitude setpoint (p.u.)
MBASE       = 7;    %% mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS  = 8;    %% status, 1 - machine in service, 0 - machine out of service
PMAX        = 9;    %% Pmax, maximum real power output (MW)
PMIN        = 10;   %% Pmin, minimum real power output (MW)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_PMAX     = 11;   %% Kuhn-Tucker multiplier on upper Pg limit (u/MW)
MU_PMIN     = 12;   %% Kuhn-Tucker multiplier on lower Pg limit (u/MW)
MU_QMAX     = 13;   %% Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
MU_QMIN     = 14;   %% Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

return;

function [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus
%IDX_BUS   Defines constants for named column indices to bus matrix.
%   [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%   VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pd = bus(4, PD);       % get the real power demand at bus 4
%    bus(:, VMIN) = 0.95;   % set the min voltage magnitude to 0.95 at all buses
% 
%   The index, name and meaning of each column of the bus matrix is given
%   below:
%
%   columns 1-13 must be included in input matrix (in case file)
%    1  BUS_I       bus number (1 to 29997)
%    2  BUS_TYPE    bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
%    3  PD          Pd, real power demand (MW)
%    4  QD          Qd, reactive power demand (MVAr)
%    5  GS          Gs, shunt conductance (MW at V = 1.0 p.u.)
%    6  BS          Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
%    7  BUS_AREA    area number, 1-100
%    8  VM          Vm, voltage magnitude (p.u.)
%    9  VA          Va, voltage angle (degrees)
%    10 BASE_KV     baseKV, base voltage (kV)
%    11 ZONE        zone, loss zone (1-999)
%    12 VMAX        maxVm, maximum voltage magnitude (p.u.)
%    13 VMIN        minVm, minimum voltage magnitude (p.u.)
%   
%   columns 14-17 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    14 LAM_P       Lagrange multiplier on real power mismatch (u/MW)
%    15 LAM_Q       Lagrange multiplier on reactive power mismatch (u/MVAr)
%    16 MU_VMAX     Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
%    17 MU_VMIN     Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
% 
%   additional constants, used to assign/compare values in the BUS_TYPE column
%    1  PQ    PQ bus
%    2  PV    PV bus
%    3  REF   reference bus
%    4  NONE  isolated bus

%   MATPOWER
%   $Id: idx_bus.m,v 1.5 2004/12/15 22:46:28 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define bus types
PQ      = 1;
PV      = 2;
REF     = 3;
NONE    = 4;

%% define the indices
BUS_I       = 1;    %% bus number (1 to 29997)
BUS_TYPE    = 2;    %% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
PD          = 3;    %% Pd, real power demand (MW)
QD          = 4;    %% Qd, reactive power demand (MVAr)
GS          = 5;    %% Gs, shunt conductance (MW at V = 1.0 p.u.)
BS          = 6;    %% Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
BUS_AREA    = 7;    %% area number, 1-100
VM          = 8;    %% Vm, voltage magnitude (p.u.)
VA          = 9;    %% Va, voltage angle (degrees)
BASE_KV     = 10;   %% baseKV, base voltage (kV)
ZONE        = 11;   %% zone, loss zone (1-999)
VMAX        = 12;   %% maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
VMIN        = 13;   %% minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
LAM_P       = 14;   %% Lagrange multiplier on real power mismatch (u/MW)
LAM_Q       = 15;   %% Lagrange multiplier on reactive power mismatch (u/MVAr)
MU_VMAX     = 16;   %% Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
MU_VMIN     = 17;   %% Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)

return;

function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch
%IDX_BRCH   Defines constants for named column indices to branch matrix.
%   [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
%   RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    branch(4, BR_STATUS) = 0;              % take branch 4 out of service
%    Ploss = branch(:, PF) + branch(:, PT); % compute real power loss vector
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-11 must be included in input matrix (in case file)
%    1  BUS_I       bus number (1 to 29997)
%    1  F_BUS       f, from bus number
%    2  T_BUS       t, to bus number
%    3  BR_R        r, resistance (p.u.)
%    4  BR_X        x, reactance (p.u.)
%    5  BR_B        b, total line charging susceptance (p.u.)
%    6  RATE_A      rateA, MVA rating A (long term rating)
%    7  RATE_B      rateB, MVA rating B (short term rating)
%    8  RATE_C      rateC, MVA rating C (emergency rating)
%    9  TAP         ratio, transformer off nominal turns ratio
%    10 SHIFT       angle, transformer phase shift angle (degrees)
%    11 BR_STATUS   initial branch status, 1 - in service, 0 - out of service
%
%   columns 12-15 are added to matrix after power flow or OPF solution
%   they are typically not present in the input matrix
%    12 PF          real power injected at "from" bus end (MW)
%    13 QF          reactive power injected at "from" bus end (MVAr)
%    14 PT          real power injected at "to" bus end (MW)
%    15 QT          reactive power injected at "to" bus end (MVAr)
%
%   columns 14-17 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    16 MU_SF       Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
%    17 MU_ST       Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)

%   MATPOWER
%   $Id: idx_brch.m,v 1.5 2004/12/15 22:46:27 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
F_BUS       = 1;    %% f, from bus number
T_BUS       = 2;    %% t, to bus number
BR_R        = 3;    %% r, resistance (p.u.)
BR_X        = 4;    %% x, reactance (p.u.)
BR_B        = 5;    %% b, total line charging susceptance (p.u.)
RATE_A      = 6;    %% rateA, MVA rating A (long term rating)
RATE_B      = 7;    %% rateB, MVA rating B (short term rating)
RATE_C      = 8;    %% rateC, MVA rating C (emergency rating)
TAP         = 9;    %% ratio, transformer off nominal turns ratio
SHIFT       = 10;   %% angle, transformer phase shift angle (degrees)
BR_STATUS   = 11;   %% initial branch status, 1 - in service, 0 - out of service

%% included in power flow solution, not necessarily in input
PF          = 12;   %% real power injected at "from" bus end (MW)       (not in PTI format)
QF          = 13;   %% reactive power injected at "from" bus end (MVAr) (not in PTI format)
PT          = 14;   %% real power injected at "to" bus end (MW)         (not in PTI format)
QT          = 15;   %% reactive power injected at "to" bus end (MVAr)   (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_SF       = 16;   %% Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST       = 17;   %% Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)

return;
