function [Hist_True, Z, t_vec, Leak_True, P_GTU] = dse_3_gen_data_leak(Nodes, Pipes, Compressors, Sys, Leaks, GTU)
% DSE_3_GEN_DATA_LEAK - Generate synthetic measurement data with leaks
% Simulates gas network dynamics and produces noisy measurements
% Outputs:
%   Hist_True  - True state history [N x nS] (density + flow)
%   Z          - Noisy measurements [N x nS]
%   t_vec      - Time vector [hours]
%   Leak_True  - Actual leak rates per pipe [N x nP]
%   P_GTU      - GTU power output [N x nG]
LHV = 50;    % Lower heating value of gas [MJ/kg]

% Initialize dimensions and storage arrays
N = floor(Sys.Hours * 3600 / Sys.dt);
nN = height(Nodes);
nP = height(Pipes);
nS = nN + nP;
nG = height(GTU);

Hist_True = zeros(N,nS);
Z = zeros(N,nS);
Leak_True = zeros(N,nP);
P_GTU = zeros(N,max(nG,1));
t_vec = (1:N) * Sys.dt / 3600;

x = dse_steady_solver(Nodes, Pipes, Compressors, Sys);

% Time-stepping Simulation Loop
for t = 1:N
    hr = t_vec(t);
    sec = t * Sys.dt;
    
    % --- Compute time-varying load profile
    lf = 1.0 + 0.05*sin(2*pi*(hr-8)/24);
    load_now = Nodes.BaseLoad .* lf;
    
    % --- Compute GTU power and gas consumption ---
    elf = 0.7 + 0.3*(sin(2*pi*(hr-6)/12).^2 .* (hr>=6 & hr<=22));  
    
    for g = 1:nG
        if GTU.Status(g) ~= 1, continue; end
        cap = GTU.Capacity_MW(g);
        pmin = GTU.MinLoad(g)/100;
        pmax = GTU.MaxLoad(g)/100;
        pwr = cap * (pmin + (pmax-pmin)*elf);
        P_GTU(t,g) = pwr;
        gas = pwr * GTU.HeatRate(g) * 1000 / 3600 / LHV;
        load_now(GTU.NodeID(g)) = gas;
    end
    
    % --- Compute leak-induced loads ---
    leak_load = zeros(nN,1);
    leak_pipe = zeros(nP,1);
    
    for lk = 1:height(Leaks)
        pid = Leaks.PipeID(lk);
        if sec >= Leaks.StartTime_s(lk) && sec <= Leaks.EndTime_s(lk)
            rate = Leaks.LeakRate_kg_s(lk);
            pos = Leaks.Position(lk);
            leak_pipe(pid) = leak_pipe(pid) + rate;
            
             n1 = Pipes.From(pid);
            n2 = Pipes.To(pid);
            if Nodes.Type(n1) ~= 1
                leak_load(n1) = leak_load(n1) + rate*(1-pos);
            end
            if Nodes.Type(n2) ~= 1
                leak_load(n2) = leak_load(n2) + rate*pos;
            end
        end
    end
    Leak_True(t,:) = leak_pipe';
    
    total_load = load_now + leak_load;
    

    % --- Create state noise ---
    noise_state_p = sqrt(1.5e-4) * randn(1,nN) * 1e5 / Sys.c2;
    noise_state_m = sqrt(4e-3) * randn(1,nP);

    [A, B, u] = dse_2_build_model(x, Nodes, Pipes, Compressors, Sys, total_load);
    x = A \ (B*x + u + [noise_state_p'; noise_state_m']);
    Hist_True(t,:) = x';
    
    % --- Add measurement noise ---
    noise_p = sqrt(0.04) * randn(1,nN) * 1e5 / Sys.c2;
    noise_m = sqrt(2) * randn(1,nP);
    Z(t,:) = Hist_True(t,:) + [noise_p, noise_m];
end
end