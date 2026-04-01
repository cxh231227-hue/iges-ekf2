function H_E = power_build_HE(PS)
% POWER_BUILD_HE  Build PMU measurement matrix H_E per Chen (2022) Appendix A.
%
% The power system state vector is:
%   x_E = [e1, f1, e2, f2, ..., e39, f39]'   (2*nB x 1)
% where e_i = Re(V_i), f_i = Im(V_i) in rectangular coordinates.
%
% The measurement vector z_E has three blocks:
%   z_V  : voltage phasors at PMU buses         (2*nZB rows)
%   z_IB : branch current phasors at PMU lines  (2*nZC rows)
%   z_IN : bus injected current phasors at PMU  (2*nZB rows)
%   z_E  = [z_V; z_IB; z_IN],  total nZE = 4*nZB + 2*nZC rows
%
% H_E = [H_V; H_IB; H_IN]   size: nZE x 2*nB
%
% Branch model: standard pi-model including off-nominal tap transformers.
% All angle shifts in case39 are zero, so t = real scalar tau.
%
% Input:
%   PS  - struct from power_init_case39
% Output:
%   H_E - measurement matrix (nZE x 2*nB), real-valued

nB          = PS.nB;
G           = PS.G;
B           = PS.B;
pmu_buses   = PS.pmu_buses;
pmu_branches= PS.pmu_branches;
branch      = PS.branch;
nZB         = PS.nZB;   % 21
nZC         = PS.nZC;   % 29
nZE         = PS.nZE;   % 142

nE = 2 * nB;   % 78

H_E = zeros(nZE, nE);

% ================================================================
% Block 1: H_V  (rows 1 : 2*nZB)
% Voltage phasor measurement at PMU bus b_M = pmu_buses(m):
%   z_V(2m-1) = e_i  ->  H_V(2m-1, 2i-1) = 1
%   z_V(2m)   = f_i  ->  H_V(2m,   2i)   = 1
% Eq. (38) in Appendix A.
% ================================================================
for m = 1:nZB
    bus_i = pmu_buses(m);
    H_E(2*m-1, 2*bus_i-1) = 1;   % e_i coefficient
    H_E(2*m,   2*bus_i)   = 1;   % f_i coefficient
end

row_offset_IB = 2*nZB;   % H_IB starts here

% ================================================================
% Block 2: H_IB  (rows 2*nZB+1 : 2*nZB+2*nZC)
% Branch current phasor measurement for PMU branch [fb, tb].
%
% For a branch (fbus=i, tbus=j) with off-nominal tap ratio tau
% (Matpower: ratio=0 means tau=1), the pi-model gives:
%
%   Y_ff = (y_s + j*b_c/2) / tau^2
%   Y_ft = -y_s / tau
%   Y_tf = -y_s / tau           (angle shift = 0)
%   Y_tt =  y_s + j*b_c/2
%
% Current at from-end (bus i):  I_ij = Y_ff*V_i + Y_ft*V_j
% Current at to-end   (bus j):  I_ji = Y_tf*V_i + Y_tt*V_j
%
% H_IB row for measurement of I flowing from bus 'a' to bus 'b':
%   If a=fbus, b=tbus: use Y_ff (row_a) and Y_ft (row_b)
%   If a=tbus, b=fbus: use Y_tf (row_a=tbus, col_fbus) and Y_tt (row_b=tbus, col_tbus)
%     -> I at tbus = Y_tf*V_fbus + Y_tt*V_tbus
%     -> Real: Re(Y_tf)*e_fbus - Im(Y_tf)*f_fbus + Re(Y_tt)*e_tbus - Im(Y_tt)*f_tbus
%
% For complex Y acting on V = e + j*f:
%   Re(Y*V) ->  Re(Y)*e - Im(Y)*f
%   Im(Y*V) ->  Im(Y)*e + Re(Y)*f
% ================================================================
for m = 1:nZC
    fb = pmu_branches(m, 1);
    tb = pmu_branches(m, 2);
    row_real = row_offset_IB + 2*m - 1;
    row_imag = row_offset_IB + 2*m;

    % Find matching branch in Matpower branch table
    br_idx = find( (branch(:,1)==fb & branch(:,2)==tb) | ...
                   (branch(:,1)==tb & branch(:,2)==fb) );

    if isempty(br_idx)
        warning('power_build_HE: no branch found for PMU branch [%d,%d]', fb, tb);
        continue;
    end
    br_idx = br_idx(1);

    fbus = branch(br_idx, 1);
    tbus = branch(br_idx, 2);
    r    = branch(br_idx, 3);
    x    = branch(br_idx, 4);
    b_c  = branch(br_idx, 5);   % total charging susceptance
    tau  = branch(br_idx, 9);
    if tau == 0, tau = 1; end   % Matpower: 0 means nominal (1.0)

    % Series admittance
    y_s = 1 / (r + 1j*x);

    % Pi-model admittances (angle shift = 0 for case39)
    Y_ff = (y_s + 1j*b_c/2) / tau^2;
    Y_ft = -y_s / tau;
    Y_tf = -y_s / tau;
    Y_tt =  y_s + 1j*b_c/2;

    if fbus == fb
        % Measurement is at from-end: I = Y_ff*V_fb + Y_ft*V_tb
        Y_a = Y_ff;   bus_a = fb;
        Y_b = Y_ft;   bus_b = tb;
    else
        % Measurement is at to-end: I = Y_tf*V_fbus + Y_tt*V_tbus
        % but measured in direction [fb=tbus -> tb=fbus]
        Y_a = Y_tt;   bus_a = fb;   % tbus side
        Y_b = Y_tf;   bus_b = tb;   % fbus side
    end

    % Fill real part row
    H_E(row_real, 2*bus_a-1) =  real(Y_a);
    H_E(row_real, 2*bus_a)   = -imag(Y_a);
    H_E(row_real, 2*bus_b-1) =  real(Y_b);
    H_E(row_real, 2*bus_b)   = -imag(Y_b);

    % Fill imag part row
    H_E(row_imag, 2*bus_a-1) =  imag(Y_a);
    H_E(row_imag, 2*bus_a)   =  real(Y_a);
    H_E(row_imag, 2*bus_b-1) =  imag(Y_b);
    H_E(row_imag, 2*bus_b)   =  real(Y_b);
end

row_offset_IN = 2*nZB + 2*nZC;   % H_IN starts here

% ================================================================
% Block 3: H_IN  (rows 2*nZB+2*nZC+1 : nZE)
% Bus injected current phasor at PMU bus b_M = pmu_buses(m).
%
% I_inj = Ybus * V   (full bus injection)
% I_b_real = sum_j [ G_bj*e_j - B_bj*f_j ]
% I_b_imag = sum_j [ B_bj*e_j + G_bj*f_j ]
%
% H_IN(2m-1, 2j-1) =  G_ij  (real injection, e_j coefficient)
% H_IN(2m-1, 2j)   = -B_ij  (real injection, f_j coefficient)
% H_IN(2m,   2j-1) =  B_ij  (imag injection, e_j coefficient)
% H_IN(2m,   2j)   =  G_ij  (imag injection, f_j coefficient)
% Eq. (40) in Appendix A.
% ================================================================
for m = 1:nZB
    bus_i    = pmu_buses(m);
    row_real = row_offset_IN + 2*m - 1;
    row_imag = row_offset_IN + 2*m;

    for j = 1:nB
        H_E(row_real, 2*j-1) =  G(bus_i, j);
        H_E(row_real, 2*j)   = -B(bus_i, j);
        H_E(row_imag, 2*j-1) =  B(bus_i, j);
        H_E(row_imag, 2*j)   =  G(bus_i, j);
    end
end

fprintf('power_build_HE: H_E size = %d x %d\n', size(H_E,1), size(H_E,2));
fprintf('  H_V  rows 1..%d       (voltage phasors at %d PMU buses)\n', 2*nZB, nZB);
fprintf('  H_IB rows %d..%d   (branch currents at %d PMU branches)\n', ...
        2*nZB+1, 2*nZB+2*nZC, nZC);
fprintf('  H_IN rows %d..%d  (injected currents at %d PMU buses)\n', ...
        2*nZB+2*nZC+1, nZE, nZB);
end
