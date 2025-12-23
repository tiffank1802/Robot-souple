% TP robot souple - controle non-causal bonus
close all;
clear all;

nstep = 1000; % Total time steps
delta = .01; % Time step
mu    = 1e-1; % Regularization
m = 100; % Prediction horizon
n = 10; % Update every n steps

beta = .25;
gamma = .5;

L = 500; a=5; b=5; E=70000; rho=2700e-9; cy=1e-4; nx=10; nA=floor(nx/2);
S=a*b; Ix=a*b^3/12; dx=L/nx;
[Kfull,Cfull,Mfull] = neb_beam_matrices(nx,dx,E,Ix,rho,S,cy);
K = Kfull; K([1,2],:)=[]; K(:,[1,2])=[];
C = Cfull; C([1,2],:)=[]; C(:,[1,2])=[];
M = Mfull; M([1,2],:)=[]; M(:,[1,2])=[];
ndo1 = size(M,1);
induA = 2*nA-1; induB = 2*nx-1;

% Reference
t0 = 1.0;
time_full = 0:delta:delta*nstep;
u_ref = double(time_full >= t0);

% Initial state
u0 = zeros(ndo1,1); v0 = zeros(ndo1,1); a0 = M \ (-C*v0 - K*u0);

% Storage
u_nc = zeros(ndo1, nstep+1);
v_nc = zeros(ndo1, nstep+1);
a_nc = zeros(ndo1, nstep+1);
f_nc = zeros(ndo1, nstep+1);
u_nc(:,1) = u0; v_nc(:,1) = v0; a_nc(:,1) = a0;

current_step = 1;
while current_step <= nstep
    % Horizon for this window
    horizon_end = min(current_step + m - 1, nstep);
    n_horizon = horizon_end - current_step + 1;
    
    % Desired displacement for horizon
    co = zeros(ndo1, n_horizon);
    ref_horizon = u_ref(current_step:horizon_end);
    co(induB, :) = ref_horizon;
    
    % Use current state as initial (assume perfect or from Kalman)
    u_init = u_nc(:, current_step);
    v_init = v_nc(:, current_step);
    a_init = a_nc(:, current_step);
    
    % Compute optimal command for horizon
    nopt = 50; stagEps = 1e-12;  % Very small to ensure computation even if initial residual is small
    [fop_h, ress] = forthodir(M, C, K, co, u_init, v_init, n_horizon, nopt, stagEps, delta, beta, gamma, mu, false, induB, induA);
    
    % Apply first min(n, n_horizon) steps
    steps_to_apply = min(n, n_horizon);
    for k = 1:steps_to_apply
        step_idx = current_step + k - 1;
        if step_idx > nstep, break; end
        f_nc(induA, step_idx) = fop_h(induA, k);
        % Simulate one step
        [u_nc(:,step_idx+1), v_nc(:,step_idx+1), a_nc(:,step_idx+1)] = newmark1stepMRHS(...
            M, C, K, f_nc(:,step_idx), u_nc(:,step_idx), v_nc(:,step_idx), a_nc(:,step_idx), delta, beta, gamma);
    end
    
    current_step = current_step + steps_to_apply;
end

% Plot
figure;
plot(time_full, u_ref, 'k--');
hold on;
plot(time_full, u_nc(induB,:), 'b');
legend('Reference', 'Non-causal control');
title('Non-causal control tracking');
xlabel('Time (s)'); ylabel('Tip displacement (m)');