function [Hist_Est, Load_Used, R_History, Dtt_History] = dse_chen_ekf(Z, Nodes, Pipes, Compressors, Sys, t_vec, GTU,mismatch)
% DSE_CHEN_EKF - Chen robust Kalman filter with innovation-based scaling
%
% Inputs:
%   Z           - Measurement matrix [N x nS]
%   Nodes       - Node table
%   Pipes       - Pipe table
%   Compressors - Compressor table
%   Sys         - System parameters
%   t_vec       - Time vector [hours]
%   GTU         - Gas turbine unit data (optional)
%
% Outputs:
%   Hist_Est    - Estimated state history [N x nS]
%   Load_Used   - Load profile used [N x nN]
%   Mu_History  - [NEW] mu_adaptive diagonal per step [N x nS]

LHV = 50;

% EKF parameters
Q_p = 0.15;     % Process noise for pressure
Q_m = 0.4;      % Process noise for flow
R_p = 1;      % Measurement noise for pressure
R_m = 10;       % Measurement noise for flow

%adj_rate=.5;
%Q_p=mismatch*Q_p;
if(mismatch<1)
    R_m=Q_m;
%R_p=adj_rate*R_p;
    R_p=Q_p;
elseif(mismatch>1)
    Q_m=R_m;
    Q_p=R_p;
end
%R_p=adj_rate*R_p;
%R_m=adj_rate*R_m;

%adj_rate=1.5;
%Q_p=adj_rate*Q_p;
%Q_m=adj_rate*Q_m;
%R_p=adj_rate*R_p;
%R_m=adj_rate*R_m;

% Chen robust filter parameters
m_w = 25;           % Sliding window length
mu_min = 1;       % Lower bound for mu
mu_max = 1000;     % Upper bound for mu
regularization = 1e-6;

% Initialize dimensions
N = length(t_vec);
nN = height(Nodes);
nP = height(Pipes);
nS = nN + nP;
nZ = nS;

has_gtu = nargin >= 7 && ~isempty(GTU) && height(GTU) > 0;

% Initialize state estimate
init_n = min(5, N);
x_k = mean(Z(1:init_n,:), 1)';

% Initialize covariance
P_k = eye(nS);
P_k(1:nN, 1:nN) = 10 * eye(nN);
P_k(nN+1:end, nN+1:end) = 25 * eye(nP);

% Process and measurement noise covariance
Q = eye(nS);
Q(1:nN, 1:nN) = Q_p * eye(nN);
Q(nN+1:end, nN+1:end) = Q_m * eye(nP);

R = eye(nS);
R(1:nN, 1:nN) = R_p * eye(nN);
R(nN+1:end, nN+1:end) = R_m * eye(nP);

% Innovation buffer (circular)
innovation_buffer = zeros(m_w, nZ);
buffer_idx = 1;
buffer_count = 0;

% Output storage
Hist_Est = zeros(N, nS);
Load_Used = zeros(N, nN);
R_History   = zeros(N, nS);
Dtt_History = zeros(N, nS);

% Main loop
for k = 1:N
    z_k = Z(k, :)';
    hr = t_vec(k);
    lf = 1.0 + 0.05 * sin(2*pi*(hr-8)/24);
    load_k = Nodes.BaseLoad * lf;
    
    if has_gtu
        elf = 0.7 + 0.3 * (sin(2*pi*(hr-6)/12).^2 * (hr>=6 && hr<=22));
        for g = 1:height(GTU)
            if GTU.Status(g) ~= 1, continue; end
            cap = GTU.Capacity_MW(g);
            pmin = GTU.MinLoad(g) / 100;
            pmax = GTU.MaxLoad(g) / 100;
            pwr = cap * (pmin + (pmax - pmin) * elf);
            gas = pwr * GTU.HeatRate(g) * 1000 / 3600 / LHV;
            load_k(GTU.NodeID(g)) = gas;
        end
    end
    
    % Prediction step
    [A, B, u] = dse_2_build_model(x_k, Nodes, Pipes, Compressors, Sys, load_k);
    x_pred = A \ (B * x_k + u);
    F = A \ B;
    P_pred = F * P_k * F' + Q;
    P_pred = (P_pred + P_pred') / 2;
   
    % Compute innovation
    H = eye(nS);
    e_k = z_k - H * x_pred;
    
    % Update innovation buffer
    innovation_buffer(buffer_idx, :) = e_k';
    buffer_idx = mod(buffer_idx, m_w) + 1;
    buffer_count = min(buffer_count + 1, m_w);
   
    % Innovation covariance adaptation
    if buffer_count >= m_w
        % Compute empirical innovation covariance
        P_F_empirical = zeros(nZ, nZ);
        for i = 1:m_w
            e_i = innovation_buffer(i, :)';
            P_F_empirical = P_F_empirical + e_i * e_i';
        end
        P_F_empirical = P_F_empirical / m_w;
        P_F_empirical = (P_F_empirical + P_F_empirical') / 2;
        P_F_empirical = P_F_empirical + regularization * eye(nZ);
      
        % Theoretical innovation covariance
        P_F_theoretical = H * P_pred * H';
        P_F_theoretical = (P_F_theoretical + P_F_theoretical') / 2;

        % Compute mu values
        diff_P = P_F_empirical - P_F_theoretical;
        mu_diag_raw = zeros(nZ, 1);
        
        for i = 1:nS
            mu_diag_raw(i) = diff_P(i,i) / R(i,i);
        end
        
        % Handle invalid values
        mu_diag_raw(isnan(mu_diag_raw)) = 1;
        mu_diag_raw(isinf(mu_diag_raw)) = 10000;
        
        % Apply constraints
        mu_diag = max(mu_min, min(mu_max, mu_diag_raw));
        mu_adaptive = diag(mu_diag);
    else
        % Warmup phase: use standard KF
        mu_adaptive = eye(nZ);
    end

    % [NEW] Record mu_adaptive diagonal for this step
    % (kept internally for covariance scaling, not returned)

    % Filter step with adaptive R
    S = H * P_pred * H' + mu_adaptive * R;
    S = (S + S') / 2;
    R_History(k, :)   = diag(R)';
    Dtt_History(k, :) = diag(S)';
    K = P_pred * H' / S;
    x_k = x_pred + K * e_k;
    
    % Joseph form covariance update
    %IKH = eye(nS) - K * H;
    %P_k = IKH * P_pred * IKH' + K * (mu_adaptive * R) * K';
    P_k = P_pred - K * H * P_pred;
    P_k = (P_k + P_k') / 2 + 1e-4 * eye(nS);
    
    % Physical constraints
    x_k(1:nN) = max(x_k(1:nN), 0.5);
    
    % Store results
    Hist_Est(k, :) = x_k';
    Load_Used(k, :) = load_k';
end

% Post-processing
for i = 1:nP
    Hist_Est(:, nN+i) = movmean(Hist_Est(:, nN+i), 3);
end

end