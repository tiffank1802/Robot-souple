% TP robot souple - comparaison optimal vs PID avec/sans perturbation
close all;
clear all;

nstep = 1000; % Nb time steps
delta = .01; % Time step
q     = .0001; % Error on model
r     = .0001; % Error on measurement

mu    = 1e-1; % Optimal control regularity parameter

Kp    = 5e-1; % PID parameters
Ki    = 1e-1;
Kd    = 1e-2;

beta = .25;
gamma = .5; % Newmark parameters

L     = 500; % Beam length
a     = 5;
b     = 5; % Beam section dimensions
E     = 70000; % Young modulus
rho   = 2700e-9; % Mass
cy    = 1e-4; % Damping coefficient
nx    = 10; % total nb of beam elements
nA    = floor(nx/2); % nb of beam elements to point A

S  = a*b;   % Section
Ix = a*b^3/12; % Moment

dx = L/nx;
time = delta:delta:delta*nstep;
time0 = [0,time];

isMatlab = true;
vseed1 = 0; vseed2 = 1;

% Get indices of points
npt = nx + 1; ndof = 2*npt; % nb dofs before BCs
ndo1 = ndof-2; % nb of dofs after BCs
induA = 2*nA-1;
induB = 2*nx-1;

% Perturbations
fpert = generateNoiseTemporal( time0, 1, q, vseed1, isMatlab ); % Perturbation force

% Operators and BCs
[Kfull,Cfull,Mfull] = neb_beam_matrices( nx, dx, E, Ix, rho, S, cy ); % Stiffness & Co.
K = Kfull; K(:,[1,2]) = []; K([1,2],:) = [];
C = Cfull; C(:,[1,2]) = []; C([1,2],:) = [];
M = Mfull; M(:,[1,2]) = []; M([1,2],:) = [];

% Initial conditions
u0 = zeros(ndo1,1); v0 = zeros(ndo1,1);

% Consigne: desired displacement at tip
t0 = 1.0;
co = zeros(ndo1, nstep);
start_idx = floor(t0 / delta) + 1;
co(induB, start_idx:end) = 1;

% Parameters for forthodir
niter = nstep;
nopt = 100;
stagEps = 1e-4;
picture = false;
indices = induB;
action = induA;

% Compute optimal command without perturbation
[fop, ress] = forthodir(M, C, K, co, u0, v0, niter, nopt, stagEps, delta, beta, gamma, mu, picture, indices, action);

fprintf('Optimal command computed, iterations: %d, final res: %.2e\n', length(ress), ress(end));

% Simulation optimal without perturbation
f_sim_opt = zeros(ndo1, nstep+1);
f_sim_opt(induA, 2:end) = fop(induA, :);

u_opt_no_pert = simulate_system(M, C, K, f_sim_opt, u0, v0, delta, beta, gamma);

% Simulation optimal with perturbation
f_sim_opt_pert = f_sim_opt;
f_sim_opt_pert(induA, :) = f_sim_opt(induA, :) + fpert';

u_opt_pert = simulate_system(M, C, K, f_sim_opt_pert, u0, v0, delta, beta, gamma);

% Simulation PID without perturbation
u_pid_no_pert = simulate_pid(M, C, K, u0, v0, delta, beta, gamma, Kp, Ki, Kd, t0, induA, induB, nstep, zeros(1,nstep+1));

% Simulation PID with perturbation
u_pid_pert = simulate_pid(M, C, K, u0, v0, delta, beta, gamma, Kp, Ki, Kd, t0, induA, induB, nstep, fpert);

% Plot comparison
figure; hold on;
ref = double(time0 >= t0);
plot(time0, ref, 'k--', 'linewidth', 3);
plot(time0, u_opt_no_pert(induB,:), 'b-', 'linewidth', 2);
plot(time0, u_opt_pert(induB,:), 'b--', 'linewidth', 2);
plot(time0, u_pid_no_pert(induB,:), 'r-', 'linewidth', 2);
plot(time0, u_pid_pert(induB,:), 'r--', 'linewidth', 2);
legend('Reference', 'Optimal no pert', 'Optimal with pert', 'PID no pert', 'PID with pert');
title('Tip displacement comparison');
xlabel('Time (s)'); ylabel('Displacement (m)');

function u_full = simulate_system(M, C, K, f_sim, u0, v0, delta, beta, gamma)
    nstep = size(f_sim,2) - 1;
    u = zeros(size(u0,1), nstep+1);
    v = zeros(size(u0,1), nstep+1);
    a = zeros(size(u0,1), nstep+1);
    a0 = M \ (f_sim(:,1) - C*v0 - K*u0);
    a(:,1) = a0;
    u(:,1) = u0;
    v(:,1) = v0;
    for step = 2:nstep+1
        [u(:,step), v(:,step), a(:,step)] = newmark1stepMRHS(...
            M, C, K, f_sim(:,step), u(:,step-1), v(:,step-1), a(:,step-1), delta, beta, gamma);
    end
    u_full = u;
end

function u_full = simulate_pid(M, C, K, u0, v0, delta, beta, gamma, Kp, Ki, Kd, t0, induA, induB, nstep, fpert)
    u_ref = double([0, delta:delta:delta*nstep] >= t0);
    integ_e = 0;
    prev_e = 0;
    u = zeros(size(u0,1), nstep+1);
    v = zeros(size(u0,1), nstep+1);
    a = zeros(size(u0,1), nstep+1);
    a0 = M \ (- C*v0 - K*u0);
    a(:,1) = a0;
    u(:,1) = u0;
    v(:,1) = v0;
    for step = 2:nstep+1
        uB_meas = u(induB, step-1);
        e = u_ref(step) - uB_meas;
        integ_e = integ_e + e*delta;
        de = (e - prev_e)/delta;
        Fpid = Kp*e + Ki*integ_e + Kd*de;
        f_k = zeros(size(u0,1),1);
        f_k(induA) = Fpid + fpert(step);
        [u(:,step), v(:,step), a(:,step)] = newmark1stepMRHS(...
            M, C, K, f_k, u(:,step-1), v(:,step-1), a(:,step-1), delta, beta, gamma);
        prev_e = e;
    end
    u_full = u;
end

% Simulation with optimal force
f_sim = zeros(ndo1, nstep+1);
f_sim(induA, 2:end) = fop(induA, :); % Apply optimal force at induA from step 2 to end

u = zeros(ndo1,nstep+1);
v = zeros(ndo1,nstep+1);
a = zeros(ndo1,nstep+1);
a0 = M \ (f_sim(:,1) - C*v0 - K*u0);
a(:,1) = a0;
u(:,1) = u0;
v(:,1) = v0;

for step = 2:nstep+1
    [u(:,step), v(:,step), a(:,step)] = newmark1stepMRHS(...
        M, C, K, f_sim(:,step), u(:,step-1), v(:,step-1), a(:,step-1), delta, beta, gamma);
end

