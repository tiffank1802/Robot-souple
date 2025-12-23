% TP robot souple - etude impact de mu
close all;
clear all;

% Define mu values to test
mu_values = [0.01, 0.1, 1];

% Store results
results = struct();

nstep = 1000; % Nb time steps
delta = .01; % Time step

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

% Get indices of points
npt = nx + 1; ndof = 2*npt; % nb dofs before BCs
ndo1 = ndof-2; % nb of dofs after BCs
induA = 2*nA-1;
induB = 2*nx-1;

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
nopt = 200; % Increased
stagEps = 1e-6; % Tighter
picture = false;
indices = induB;
action = induA;

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    % Call forthodir
    [fop, ress] = forthodir(M, C, K, co, u0, v0, niter, nopt, stagEps, delta, beta, gamma, mu, picture, indices, action);
    
    % Convergence metrics
    num_iter = length(ress);
    init_res = ress(1);
    final_res = ress(end);
    res_ratio = final_res / init_res;
    
    % Simulation with optimal force
    f_sim = zeros(ndo1, nstep+1);
    f_sim(induA, 2:end) = fop(induA, :);
    
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
    
    % Quality metrics
    ref = [0, co(induB,:)];
    tracking_error = mean(abs(u(induB,:) - ref));
    max_error = max(abs(u(induB,:) - ref));
    
    % Command amplitude
    command_amplitude = max(abs(fop(induA,:)));
    
    % Store results
    results(i).mu = mu;
    results(i).num_iter = num_iter;
    results(i).res_ratio = res_ratio;
    results(i).tracking_error = tracking_error;
    results(i).max_error = max_error;
    results(i).command_amplitude = command_amplitude;
    results(i).u_tip = u(induB,:);
    results(i).ress = ress;
    
    fprintf('mu = %.2f: iter=%d, res_ratio=%.2e, track_err=%.2e, max_err=%.2e, cmd_amp=%.2e\n', ...
        mu, num_iter, res_ratio, tracking_error, max_error, command_amplitude);
end

% Results summary
fprintf('\nSummary:\n');
fprintf('mu\tIterations\tRes Ratio\tTrack Error\tMax Error\tCmd Amplitude\n');
for i = 1:length(results)
    fprintf('%.2f\t%d\t\t%.2e\t%.2e\t%.2e\t%.2e\n', ...
        results(i).mu, results(i).num_iter, results(i).res_ratio, ...
        results(i).tracking_error, results(i).max_error, results(i).command_amplitude);
end

% Plot solutions
figure; hold on;
ref = [0, co(induB,:)];
plot(time0, ref, 'k--', 'linewidth', 2);
for i = 1:length(results)
    plot(time0, results(i).u_tip, 'linewidth', 2);
end
legend(['Reference', arrayfun(@(x) sprintf('mu=%.2f', x), [results.mu], 'UniformOutput', false)]);
title('Tip displacement for different mu');
xlabel('Time (s)'); ylabel('Displacement (m)');

% Plot residuals
figure; hold on;
for i = 1:length(results)
    plot(1:length(results(i).ress), results(i).ress, 'linewidth', 2);
end
legend(arrayfun(@(x) sprintf('mu=%.2f', x), [results.mu], 'UniformOutput', false));
title('Residual convergence for different mu');
xlabel('Iteration'); ylabel('Residual norm');
set(gca, 'YScale', 'log');

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

