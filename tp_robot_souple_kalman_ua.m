% TP robot souple - Kalman only, measuring uA
close all;
clear all;

nstep = 1000; % Nb time steps
delta = .01; % Time step
q     = .0001; % Error on model
r     = .0001; % Error on measurement

beta = .25;
gamma = .5; % Newmark parameters

L     = 500; % Beam length
a     = 5;
b     = 5; % Beam section dimensions
E     = 70000; % Young modulus
rho   = 2700e-9; % Mass
cy    = 1e-4; % Damping coefficient
nx    = 10; % total nb of beam elements
nA    = floor(nx/2);%floor(9*nx/10); % nb of beam elements to point A (watch out nx has to be multiple of 2)

S  = a*b;   % Section
Ix = a*b^3/12; % Moment

dx = L/nx;
time = delta:delta:delta*nstep;
time0 = [0,time];

isMatlab = true;
vseed1 = 0; vseed2 = 1; vseed3 = 2;

% Get indices of points
npt = nx + 1; ndof = 2*npt; % nb dofs before BCs
ndo1 = ndof-2; % nb of dofs after BCs
induA = 2*nA-1;
induB = 2*nx-1;

% Perturbations
fpert = generateNoiseTemporal( time0, 1, q, vseed1, isMatlab ); % Perturbation force
mpert = generateNoiseTemporal( time0, 1, r, vseed2, isMatlab ); % Measure perturbation

% Operators and BCs
[Kfull,Cfull,Mfull] = neb_beam_matrices( nx, dx, E, Ix, rho, S, cy ); % Stiffness & Co.
K = Kfull; K(:,[1,2]) = []; K([1,2],:) = [];
C = Cfull; C(:,[1,2]) = []; C([1,2],:) = [];
M = Mfull; M(:,[1,2]) = []; M([1,2],:) = [];

% Initial conditions
u0 = zeros(ndo1,1); v0 = zeros(ndo1,1);
a0 = M \ (-C*v0 - K*u0);

% Simulation open-loop
u = zeros(ndo1,nstep+1);
v = zeros(ndo1,nstep+1);
a = zeros(ndo1,nstep+1);
a(:,1) = a0;
u(:,1) = u0;
v(:,1) = v0;

for step = 2:nstep+1
    f_k = zeros(ndo1,1);
    f_k(induA) = fpert(step); % Perturbation at A
    [u(:,step), v(:,step), a(:,step)] = newmark1stepMRHS(...
        M, C, K, f_k, u(:,step-1), v(:,step-1), a(:,step-1), delta, beta, gamma);
end

% Kalman setup
B = Bconstruct(M,C,K,nA,delta,beta,gamma);
F = Fconstruct(M,C,K,delta,beta,gamma);
H = zeros(1,3*ndo1); H(induA) = 1; % Measure uA

x_est = [u0; v0; a0]; % Initial estimate
P = zeros(3*ndo1,3*ndo1); % Initial covariance
Q = q * B * B';
R = r;

u_est = zeros(ndo1, nstep+1);
u_est(:,1) = u0;
sigma_u = zeros(3*ndo1, nstep+1);
errors_pred = [];
errors_upd = [];

for step = 2:nstep+1
    % Prediction
    f_k_pred = zeros(ndo1,1);
    f_k_pred(induA) = fpert(step); % Assume known perturbation for prediction
    x_pred = F * x_est + B * f_k_pred;
    P_pred = F * P * F' + Q;

    % Predicted uA
    uA_pred = H * x_pred;
    error_pred = u(induA, step) - uA_pred; % Effective error
    errors_pred = [errors_pred, error_pred];

    % Measurement
    y = u(induA, step) + mpert(step); % True uA + noise

    % Update
    innov = y - H * x_pred;
    S = H * P_pred * H' + R;
    K_kalman = P_pred * H' / S;
    x_est = x_pred + K_kalman * innov;
    P = (eye(3*ndo1) - K_kalman * H) * P_pred;

    % Updated uB estimate
    uB_est = x_est(induB);
    u_est(:,step) = x_est(1:ndo1);

    error_upd = u(induB, step) - uB_est; % Effective error for uB
    errors_upd = [errors_upd, error_upd];

    sigma_u(:,step) = sqrt(diag(P));
end

u_P = u_est + sigma_u(1:ndo1,:);
u_M = u_est - sigma_u(1:ndo1,:);
% Check compatibility: std of errors vs estimated std
fprintf('Prediction: Mean error %.2e, Std error %.2e, Estimated std %.2e\n', ...
    mean(abs(errors_pred)), std(errors_pred), mean(sigma_u(induA,2:end)));
fprintf('Update: Mean error %.2e, Std error %.2e, Estimated std %.2e\n', ...
    mean(abs(errors_upd)), std(errors_upd), mean(sigma_u(induB,2:end)));

% Plot
figure; hold on;
plot(time0, u(induB,:),'b','linewidth',2);
plot(time0, u_est(induB,:),'r--','linewidth',2);
plot(time0, u_P(induB,:),'k:','linewidth',1);
plot(time0, u_M(induB,:),'k:','linewidth',1);

legend('True uB', 'Estimated uB', 'Estimated uB + 1 std', 'Estimated uB - 1 std');
title('Kalman Estimation of uB from uA measurement');
xlabel('Time (s)'); ylabel('Displacement (m)');