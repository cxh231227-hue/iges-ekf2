function [Hist_Est, Load_Used] = dse_normal_ekf(Z, Nodes, Pipes, Compressors, Sys, t_vec, GTU,mismatch)
% DSE_NORMAL_EKF - Standard Extended Kalman Filter without Leak Compensation
% This is the baseline EKF without innovation feedback for comparison
% Inputs:
%   Z           - Measurement matrix [N x nS]
%   Nodes, Pipes, Compressors - Network topology
%   Sys         - System parameters
%   t_vec       - Time vector
%   GTU         - Gas Turbine Units data
% Outputs:
%   Hist_Est    - Estimated state history [N x nS]
%   Load_Used   - Load profile used (without compensation) [N x nN]

LHV = 50;

% EKF Tuning Parameters (same as adaptive EKF for fair comparison)
Q_p = 0.15;     % Process noise for pressure states
Q_m = 0.4;      % Process noise for flow states
%Q_m = 10;
R_p = 1;      % Measurement noise for pressure
R_m = 10;       % Measurement noise for flow

%Q_p = 2;     % Process noise for pressure states
%Q_m = 15;      % Process noise for flow states
%R_p = 2;      % Measurement noise for pressure
%R_m = 15;       % Measurement noise for flow


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

% 
% Initialize dimensions
N = length(t_vec);
nN = height(Nodes);
nP = height(Pipes);
nS = nN + nP;

% Parse optional arguments
has_gtu = nargin >= 7 && ~isempty(GTU) && height(GTU) > 0;


% Initialize State Estimate
init_n = min(5, N);
x_k = mean(Z(1:init_n,:), 1)';

% Initialize Covariance Matrices
P_k = eye(nS);
P_k(1:nN, 1:nN) = 10 * eye(nN);
P_k(nN+1:end, nN+1:end) = 25 * eye(nP);

Q = eye(nS);
Q(1:nN, 1:nN) = Q_p * eye(nN);
Q(nN+1:end, nN+1:end) = Q_m * eye(nP);

R = eye(nS);
R(1:nN, 1:nN) = R_p * eye(nN);
R(nN+1:end, nN+1:end) = R_m * eye(nP);

% Output storage
Hist_Est = zeros(N, nS);
Load_Used = zeros(N, nN);

% Main Standard EKF Loop (NO LEAK COMPENSATION)
for k = 1:N
    % Get measurement at time k
    z_k = Z(k, :)';
    hr = t_vec(k);
    
    % Compute Current Load Profile--
    lf = 1.0 + 0.05 * sin(2*pi*(hr-8)/24);
    load_k = Nodes.BaseLoad * lf;
    
    % Add GTU gas consumption if present
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
    
    
    % TIME UPDATE (Prediction Step)
    [A, B, u] = dse_2_build_model(x_k, Nodes, Pipes, Compressors, Sys, load_k);
    
    % Predict next state
    x_pred = A \ (B * x_k + u);
    
    % Compute state transition matrix
    F = A \ B;
    
    % Predict error covariance
    P_pred = F * P_k * F' + Q;
    
    % MEASUREMENT UPDATE (Correction Step)
    H = eye(nS);
    
    % Compute innovation
    nu_k = z_k - H * x_pred;
    
    % Innovation covariance
    S = H * P_pred * H' + R;
    
    % Kalman gain
    K = P_pred * H' / S;
    
    % State update
    x_k = x_pred + K * nu_k;
    
    % Covariance update
    P_k = (eye(nS) - K * H) * P_pred;
    P_k = (P_k + P_k') / 2 + 1e-4 * eye(nS);
  
    % Physical Constraint
    x_k(1:nN) = max(x_k(1:nN), 0.5);
    
    % Store Results
    Hist_Est(k, :) = x_k';
    Load_Used(k, :) = load_k';
end

% Post-processing: Smooth flow estimates
for i = 1:nP
    Hist_Est(:, nN+i) = movmean(Hist_Est(:, nN+i), 3);
end

end