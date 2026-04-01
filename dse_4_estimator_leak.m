function [Hist_Est, Leak_Est, Load_Compensated, R_History] = dse_4_estimator_leak(Z, Nodes, Pipes, Compressors, Sys, t_vec, GTU, Det_Signal,mismatch)
% DSE_4_ESTIMATOR_LEAK - EKF with innovation feedback for leak compensation
%
% Inputs:
%   Z           - Measurement matrix [N x nS]
%   Nodes       - Node table
%   Pipes       - Pipe table
%   Compressors - Compressor table
%   Sys         - System parameters
%   t_vec       - Time vector [hours]
%   GTU         - Gas turbine unit data (optional)
%   Det_Signal  - Detection signal from leak detector [N x 1]
%
% Outputs:
%   Hist_Est         - Estimated state history [N x nS]
%   Leak_Est         - Leak estimation at nodes [N x nN]
%   Load_Compensated - Compensated load profile [N x nN]

LHV = 50;

% EKF parameters
Q_p = 0.15;
Q_m = 0.4;
R_p = 1;
R_m = 10;

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
% Innovation feedback parameters
delta_th = 0.35;        % Dead zone [Bar]
alpha_dk = 0.95;        % Smoothing factor
kappa = 40;             % Compensation gain

% Outlier suppression parameters
k_outlier = 3;
k2_outlier = 5;
k3_outlier = 10;
c_mad = 1.483;
alpha_sigma = 0.99;
window_size = 25;

% Initialization protection period
init_protect_hours = 5.0;

% Initialize dimensions
N = length(t_vec);
nN = height(Nodes);
nP = height(Pipes);
nS = nN + nP;

has_gtu = nargin >= 7 && ~isempty(GTU) && height(GTU) > 0;
has_detector = nargin >= 8 && ~isempty(Det_Signal);

% Compute protection period steps
init_protect_steps = round(init_protect_hours * 3600 / Sys.dt);

% Initialize state
init_n = min(5, N);
x_k = mean(Z(1:init_n,:), 1)';

% Initialize covariance
P_k = eye(nS);
P_k(1:nN, 1:nN) = 10 * eye(nN);
P_k(nN+1:end, nN+1:end) = 25 * eye(nP);

Q = diag([Q_p*ones(1,nN), Q_m*ones(1,nP)]);
R = diag([R_p*ones(1,nN), R_m*ones(1,nP)]);

% Compensation state
d_k = zeros(nN, 1);

% Outlier detection state
sigma_est = ones(nS, 1);
innovation_history = zeros(window_size, nS);
init_protect_history=zeros(init_protect_hours,nS);
hist_idx = 1;
hist_count = 0;

% Output arrays
Hist_Est = zeros(N, nS);
Leak_Est = zeros(N, nN);
Load_Compensated = zeros(N, nN);
R_History = zeros(N, nS);

% Main loop
for k = 1:N
    z_k = Z(k, :)';
    hr = t_vec(k);
    
    % Compute nominal load
    lf = 1.0 + 0.05 * sin(2*pi*(hr-8)/24);
    load_k = Nodes.BaseLoad * lf;
    
    % GTU consumption
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
    
    % Apply leak compensation
    for i = 1:nN
        if Nodes.Type(i) ~= 1
            load_k(i) = load_k(i) + d_k(i);
        end
    end
    
    % EKF prediction
    [A, B, u] = dse_2_build_model(x_k, Nodes, Pipes, Compressors, Sys, load_k);
    x_pred = A \ (B * x_k + u);
    F = A \ B;
    P_pred = F * P_k * F' + Q;
    
    % Innovation
    H = eye(nS);
    nu_k = z_k - H * x_pred;
    
    % Convert to physical units
    nu_pressure = nu_k(1:nN) * Sys.c2 / 1e5;  % [Bar]
    nu_flow = nu_k(nN+1:end);                  % [kg/s]
    nu_physical = [nu_pressure; nu_flow];
    
    % Update sigma estimate
    init_protect_history(k,:)=abs(nu_physical');
    innovation_history(hist_idx, :) = abs(nu_physical');
    hist_idx = mod(hist_idx, window_size) + 1;
    hist_count = min(hist_count + 1, window_size);
    
    if k>= init_protect_steps%window_size
        med_vals = median(innovation_history, 1)';
		finite_corr=1+5/(window_size-1); % finite correction factor
        sigma_est = alpha_sigma * sigma_est + (1 - alpha_sigma) * finite_corr* c_mad * med_vals;
        %% Detect outlier
        is_outlier = abs(nu_physical) > k_outlier * sigma_est;
        is_big_outlier=abs(nu_physical) > k2_outlier * sigma_est;
        is_very_big_outlier=abs(nu_physical) > k3_outlier * sigma_est;
    elseif k >= window_size
        med_vals = median(init_protect_history(1:k,:), 1)';
		finite_corr=1+5/(k-1); % finite correction factor
        sigma_est=c_mad*finite_corr*med_vals;
        is_outlier = zeros(size(nu_physical));
        is_big_outlier=zeros(size(nu_physical));
        is_very_big_outlier=zeros(size(nu_physical));
    else
        is_outlier = zeros(size(nu_physical));
        is_big_outlier=zeros(size(nu_physical));
        is_very_big_outlier=zeros(size(nu_physical));
    end
    sigma_est(1:nN) = max(sigma_est(1:nN), 0.05);
    sigma_est(nN+1:end) = max(sigma_est(nN+1:end), 0.5);
    
    % Outlier detection
    %is_outlier = abs(nu_physical) > k_outlier * sigma_est;
    
    % Adaptive R
    R_adaptive = R;
    for m = 1:nS
        if is_outlier(m)
            R_adaptive(m, m) = abs(nu_physical(m))/k_outlier; % Huber weight function
        elseif is_big_outlier(m)
            R_adaptive(m, m) = (k3_outlier-k2_outlier)/(k_outlier*(k3_outlier/abs(nu_physical(m))-1)); %% BUG FIX: was k3_ouliter (typo in Andrew's template)
        elseif is_very_big_outlier(m)
            R_adaptive(m, m) = 1e6;
        end
    end

    
    % EKF update
    S = H * P_pred * H' + R_adaptive;
    R_History(k, :) = diag(R)';
    K = P_pred * H' / S;
    x_k = x_pred + K * nu_k;
    P_k = (eye(nS) - K * H) * P_pred;
    P_k = (P_k + P_k') / 2 + 1e-4 * eye(nS);
    
    % Check initialization period
    in_init_period = (k <= init_protect_steps);
    
    % Check detector signal
    leak_detected = false;
    if has_detector && k <= length(Det_Signal)
        leak_detected = Det_Signal(k);
    end
    
    % Update leak compensation
    for i = 1:nN
        if Nodes.Type(i) == 1  % Skip source nodes
            continue;
        end
        
        % Force zero during initialization
        if in_init_period
            d_k(i) = 0;
            continue;
        end
        
        % Freeze during outlier
        if is_outlier(i)
            d_k(i) = d_k(i) * 0.98;
            continue;
        end
        
        nu_p_i = nu_pressure(i);
        
        if leak_detected
            if nu_p_i < -delta_th
                % Negative innovation = pressure drop = leak
                leak_signal = -nu_p_i;
                d_k(i) = alpha_dk * d_k(i) + (1 - alpha_dk) * kappa * leak_signal;
                d_k(i) = max(0, min(d_k(i), 80));
            elseif nu_p_i > delta_th && d_k(i) > 0
                % Positive innovation with compensation = over-compensation
                d_k(i) = d_k(i) * 0.85;
            elseif abs(nu_p_i) <= delta_th
                % Within dead zone
                d_k(i) = d_k(i) * 0.995;
            end
        else
            % No detection: slow decay
            d_k(i) = d_k(i) * 0.95;
        end
        
        % Clear small values
        if d_k(i) < 0.1
            d_k(i) = 0;
        end
    end
    
    % Physical constraints
    x_k(1:nN) = max(x_k(1:nN), 0.5);
    
    % Store results
    Hist_Est(k, :) = x_k';
    Leak_Est(k, :) = d_k';
    Load_Compensated(k, :) = load_k';
end

% Post-processing
for i = 1:nP
    Hist_Est(:, nN+i) = movmean(Hist_Est(:, nN+i), 3);
end

end