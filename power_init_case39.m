function PS = power_init_case39()
% POWER_INIT_CASE39  Initialize IEEE 39-bus power system data for IGES coupling.
%
% Loads case39, runs power flow via Matpower, extracts admittance matrix,
% defines GTU coupling and PMU measurement configuration per Chen (2022).
%
% GTU mapping (power bus -> Nodes table row in dse_load_gaslib40 output):
%   bus 31 -> Row 31 (sink_28, Chen node 31)
%   bus 32 -> Row 24 (sink_21, Chen node 24)
%   bus 34 -> Row  4 (sink_1,  Chen node  4)
%   bus 33 -> Row 37 (innode_5,Chen node 34)  <- NOT Row34
%   bus 36 -> Row 30 (sink_27, Chen node 30)
%
% Output PS fields:
%   G, B          - real/imag parts of Ybus (39x39)
%   Ybus          - complex admittance matrix (39x39 sparse)
%   nB            - number of buses (39)
%   x_E_init      - initial power state [e1;f1;e2;f2;...;e39;f39] (78x1)
%   gtu_bus       - power bus numbers of GTUs (1x5)
%   gtu_row       - Nodes table row indices of GTUs (1x5)  <- use this
%   eta           - energy conversion coefficient [MW.s/kg]
%   alpha_S       - Holt level smoothing parameter
%   beta_S        - Holt trend smoothing parameter
%   pmu_buses     - bus indices with PMU voltage+injection measurement
%   pmu_branches  - branch endpoint pairs with current measurement (29x2)
%   branch        - branch table from case39
%   baseMVA       - system MVA base

% --- Matpower auto-detect ---
if ~exist('makeYbus', 'file')
    candidates = {
        'C:\matpower8.0', 'C:\matpower7.1', 'C:\matpower',
        fullfile(getenv('USERPROFILE'), 'matpower8.0'),
        fullfile(getenv('USERPROFILE'), 'matpower'),
        '/usr/local/matpower', '/home/matpower',
        fullfile(getenv('HOME'), 'matpower')
    };
    found = false;
    for i = 1:length(candidates)
        p = candidates{i};
        if exist(p, 'dir')
            addpath(genpath(p));
            if exist('makeYbus', 'file')
                fprintf('Matpower found at: %s\n', p);
                found = true;
                break;
            end
        end
    end
    if ~found
        error(['Matpower not found. Please run: addpath(genpath(''<your_matpower_path>''))' ...
               ' before calling power_init_case39.']);
    end
end

% --- Load case and run power flow ---
mpc    = case39;
mpopt  = mpoption('verbose', 0, 'out.all', 0);
result = runpf(mpc, mpopt);

if result.success ~= 1
    error('Matpower power flow did not converge for case39.');
end

nB = size(mpc.bus, 1);   % 39

% --- Admittance matrix ---
[Ybus, ~, ~] = makeYbus(mpc);
PS.Ybus = Ybus;
PS.G    = real(full(Ybus));   % 39 x 39
PS.B    = imag(full(Ybus));   % 39 x 39
PS.nB   = nB;

% --- Initial power state from solved power flow ---
% State vector ordering: [e1, f1, e2, f2, ..., e39, f39]  (78 x 1)
V_mag  = result.bus(:, 8);          % voltage magnitude (p.u.)
V_ang  = result.bus(:, 9) * pi/180; % voltage angle (rad)
e_init = V_mag .* cos(V_ang);       % real part
f_init = V_mag .* sin(V_ang);       % imaginary part

ef_matrix    = [e_init, f_init];    % 39 x 2
PS.x_E_init  = reshape(ef_matrix', [], 1);  % 78 x 1: [e1;f1;e2;f2;...]

PS.e_init = e_init;
PS.f_init = f_init;

% --- GTU coupling: bus -> gas node row index in Nodes table ---
% Chen (2022) Section V-A, five GTUs.
%
% IMPORTANT: Chen's official node numbers and Nodes table row indices
% differ for one GTU. Confirmed by pipe topology matching:
%
%   Power bus | Chen gas node | Nodes table Row | XML ID
%      31     |      31       |      31         | sink_28
%      32     |      24       |      24         | sink_21
%      34     |       4       |       4         | sink_1
%      33     |      34       |      37         | innode_5  <- Chen34=Row37
%      36     |      30       |      30         | sink_27
%
% Chen gas 34 = innode_5 (Row 37), connected via pipes to:
%   innode_4(Row36=ChenNode3, gas source) and sink_27(Row30=ChenNode30)
% Row 34 (innode_2) is a compressor outlet deleted in Chen's model.
PS.gtu_bus  = [31, 32, 34, 33, 36];   % power system bus numbers
PS.gtu_row  = [31, 24,  4, 37, 30];   % row indices in Nodes table (USE THIS in code)
PS.nGTU     = 5;

% Energy conversion coefficient: 20.148 MW.s/kg at 40% efficiency
% Chen (2022) Section V-A
PS.eta = 20.148;   % [MW.s/kg]

% --- Holt's exponential smoothing parameters ---
% Chen (2022) Section V-A: alpha_S = 0.5, beta_S = 0.4
PS.alpha_S = 0.5;
PS.beta_S  = 0.4;

% --- PMU measurement configuration (Chen 2022 Section V-A) ---
% Buses equipped with PMU (voltage phasor + injected current)
PS.pmu_buses = [2,4,6,9,12,15,18,21,22,25,28,30,31,32,33,34,35,36,37,38,39];

% Branch current measurements: each row is [from_bus, to_bus]
PS.pmu_branches = [
    39,  1;
    25,  2;
    30,  2;
     6,  5;
     7,  6;
    11,  6;
    31,  6;
     9,  8;
    39,  9;
    11, 10;
    13, 10;
    32, 10;
    14, 13;
    15, 14;
    17, 16;
    18, 17;
    27, 17;
    20, 19;
    33, 19;
    34, 20;
    35, 22;
    24, 23;
    36, 23;
    26, 25;
    37, 25;
    28, 26;
    29, 26;
    38, 29
];

% --- Branch data (needed for H_E construction) ---
% Columns: [fbus tbus r x b ratio angle]
br = mpc.branch;
PS.branch = br;

% --- Base ---
PS.baseMVA  = mpc.baseMVA;
PS.mpc_base = mpc;   % keep full mpc for runpf in data generation

% --- Dimension summary ---
nZB = length(PS.pmu_buses);   % 21 PMU buses
nZC = size(PS.pmu_branches, 1); % 29 branch current measurements
% Measurement vector z_E layout:
%   [V_real(nZB); V_imag(nZB); I_branch_real(nZC); I_branch_imag(nZC);
%    I_inject_real(nZB); I_inject_imag(nZB)]
% Total: 4*nZB + 2*nZC
PS.nZB  = nZB;
PS.nZC  = nZC;
PS.nZE  = 4*nZB + 2*nZC;   % = 4*21 + 2*29 = 84 + 58 = 142

fprintf('power_init_case39: %d buses, %d branches\n', nB, size(mpc.branch,1));
fprintf('  GTU buses    : [%s]\n', num2str(PS.gtu_bus));
fprintf('  GTU row idx  : [%s]\n', num2str(PS.gtu_row));
fprintf('  PMU buses (%d), branch currents (%d)\n', nZB, nZC);
fprintf('  nE = %d,  nZE = %d\n', 2*nB, PS.nZE);
end
