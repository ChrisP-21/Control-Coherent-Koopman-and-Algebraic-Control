clc
clear
close all

% mainCartPend.m
% This file computes Koopman models of the cart pendulum through gEDMD
% (with a trajectory dataset and a uniformly random dataset), DE, and CCK.
% This file also compares these models to the original nonlinear model and 
% to each other

%%

rng(26); % for repeatability

saveModels = false; % set to true to save DE, EDMD, and CCK models
saveRMSEs = false; % set to true to save RMSE values from simulation trials

%% 
tic

% DE on cart pendulum

% time step size [s]
dt = 0.02;

% cart pendulum parameters
M = 1.1; % mass of cart [kg]
m = 0.9; % mass of pendulum [kg]
g = 9.81; % gravity acceleration [m/s^2]
L = 1; % length of pendulum [m]
cd = 1; % pendulum damping coefficient [N*s/m]

% helper variables
k1 = g/L;
k2 = cd/(m*L^2);
k3 = L*(M + m);
k4 = m*L;
k5 = m/(M + m);

% defining symbolic variables for state and input
syms x1 x2 u real
x = [x1; x2];

% cart pendulum state equation
f_cont = [x2;
          (-k1*sin(x1) - k2*x2 - cos(x1)/k3 .* (u + k4*sin(x1).*x2.^2))./(1 - k5*cos(x1).^2)];

% output equation
h = [x1; x2];

% constructing Hermite polynomial basis for DE and EDMD
rho = 6; % maximum total degree
psi = hermiteBasis(x, rho, [rho, rho]);

% integration parameters
bounds = [-inf(length([x; u]), 1), inf(length([x; u]), 1)];
w = 1;
int_params = {bounds, w, "Gauss-Hermite", 20};

% compute (A, B, C) matrices for DE model
[A_DE, B_DE, C_DE] = DirectEncoding(f_cont, h, x, psi, u, "MDI-TP", int_params);

if (saveModels)
    save K_DE.mat A_DE B_DE C_DE f_cont h psi x u
end

%% 

rng(randi(50)); % for repeatability

NLS = {f_cont, h};
KSS = {A_DE, B_DE, C_DE};

TF = 30; % simulation end time [s]
Ntrials = 1000; % number of simulation trials

% randomized initial conditions for testing
x0_test = [-pi/2; -1] + [pi; 2] .* rand(length(x), Ntrials);

% randomized input waveforms for testing
u_test = cell(1, Ntrials);
for i = 1:Ntrials
    t_samp = linspace(0, TF, 16).';
    u_samp = -10 + 20*rand(16, length(u));
    u_test{i} = @(t) interp1(t_samp, u_samp, t).';
end

% plot settings
pltset = cell(1, 2);
pltset{1} = 1;
pltset{2} = {"Angular Position (CT-DE)", "Angular Velocity (CT-DE)"};

% comparing DE model to nonlinear model
RMSEs_DE = NLS_vs_KoopmanCont(NLS, KSS, psi, x, u, dt, TF, x0_test, u_test, pltset);

%%
rng(15); % for repeatability

% gEDMD with uniformly random dataset
[X, U, X_dot, Y] = collectDataCT_Uni(f_cont, h, x, u, 40e3, [-pi/2; -5], [pi/2; 5], -10, 10);
[A_EDMD_uni, B_EDMD_uni, C_EDMD_uni] = CT_EDMD(X, U, X_dot, Y, psi, x);

if (saveModels)
    save K_EDMD_uni.mat A_EDMD_uni B_EDMD_uni C_EDMD_uni f_cont h psi x u
end

%% 

% comparing gEDMD model (with uniformly random data) to nonlinear model
KSS = {A_EDMD_uni, B_EDMD_uni, C_EDMD_uni};
pltset{2} = {"Angular Position (gEDMD, uniformly random dataset)", "Angular Velocity (gEDMD, uniformly random dataset)"};
RMSEs_EDMD_uni = NLS_vs_KoopmanCont(NLS, KSS, psi, x, u, dt, TF, x0_test, u_test, pltset);

%% 
rng(34); % for repeatability

% gEDMD with trajectory dataset
[X, U, X_dot, Y] = collectDataCT_Traj(f_cont, h, x, u, 20, 200, dt, [-pi/2; -1], [pi/2; 1], -10, 10);
[A_EDMD_traj, B_EDMD_traj, C_EDMD_traj] = CT_EDMD(X, U, X_dot, Y, psi, x);

if (saveModels)
    save K_EDMD_traj.mat A_EDMD_traj B_EDMD_traj C_EDMD_traj f_cont h psi x u
end

%% 

% comparing gEDMD model to nonlinear model
KSS = {A_EDMD_traj, B_EDMD_traj, C_EDMD_traj};
pltset{2} = {"Angular Position (gEDMD, trajectory dataset)", "Angular Velocity (gEDMD, trajectory dataset)"};
RMSEs_EDMD_traj = NLS_vs_KoopmanCont(NLS, KSS, psi, x, u, dt, TF, x0_test, u_test, pltset);

%% 

% CCK on cart pendulum

% actuator states
syms p1 p2 real
p = [p1; p2];

% actuator-augmented cart pendulum state equation
b = -100; % low-pass filter pole
f_cont2 = [p2;
           -b^2*p1 + 2*b*p2 + b^2*u;
           x2;
           (-k1*sin(x1) - k2*x2 - cos(x1)/k3 .* (p1 + k4*sin(x1).*x2.^2))./(1 - k5*cos(x1).^2)];

% output equation
h = [x1; x2];

% constructing Hermite polynomial basis for CCK
rho = 6;
psi_CCK = hermiteBasis([p; x], 6, [1, 1, rho, rho]);

% integration parameters
bounds = [-inf(length([p; x; u]), 1), inf(length([p; x; u]), 1)];
w = 1;
int_params = {bounds, w, "Gauss-Hermite", 20};

% compute (A, B, C) matrices for CCK model
[A_CCK, B_CCK, C_CCK] = DirectEncoding(f_cont2, h, [p; x], psi_CCK, u, "MDI-TP", int_params);

if (saveModels)
    save K_CCK.mat A_CCK B_CCK C_CCK f_cont2 h psi_CCK p x u
end

%% 

% comparing Koopman model to nonlinear model
NLS = {f_cont2, h};
KSS = {A_CCK, B_CCK, C_CCK};
x0_test2 = [zeros(length(p), Ntrials); x0_test];
pltset{2} = {"Angular Position (CT-CCK)", "Angular Velocity (CT-CCK)"};
RMSEs_CCK = NLS_vs_KoopmanCont(NLS, KSS, psi_CCK, [p; x], u, dt, TF, x0_test2, u_test, pltset);

%%

% comparing RMSEs of different Koopman identification methods

not_nan_traj = ~isnan(RMSEs_EDMD_traj);
not_nan_uni = ~isnan(RMSEs_EDMD_uni);
not_nan_DE = ~isnan(RMSEs_DE);
not_nan_CCK = ~isnan(RMSEs_CCK);

% averages (for average prediction accuracy)
fprintf("Average RMSEs\n")
fprintf("EDMD (Trajectory): %.5f\n", mean(RMSEs_EDMD_traj(not_nan_traj)))
fprintf("EDMD (Uniform): %.5f\n", mean(RMSEs_EDMD_uni(not_nan_uni)))
fprintf("DE: %.5f\n", mean(RMSEs_DE(not_nan_DE)))
fprintf("CCK: %.5f\n\n", mean(RMSEs_CCK(not_nan_CCK)))

% standard deviations (for the consistency of prediction accuracy)
fprintf("Standard Deviation RMSEs\n")
fprintf("EDMD (Trajectory): %.5f\n", std(RMSEs_EDMD_traj(not_nan_traj)))
fprintf("EDMD (Uniform): %.5f\n", std(RMSEs_EDMD_uni(not_nan_uni)))
fprintf("DE: %.5f\n", std(RMSEs_DE(not_nan_DE)))
fprintf("CCK: %.5f\n\n", std(RMSEs_CCK(not_nan_CCK)))

% minima (for best-case-scenario accuracy)
fprintf("Minimum RMSEs\n")
fprintf("EDMD (Trajectory): %.5f\n", min(RMSEs_EDMD_traj(not_nan_traj)))
fprintf("EDMD (Uniform): %.5f\n", min(RMSEs_EDMD_uni(not_nan_uni)))
fprintf("DE: %.5f\n", min(RMSEs_DE(not_nan_DE)))
fprintf("CCK: %.5f\n\n", min(RMSEs_CCK(not_nan_CCK)))

% maxima (for worst-case-scenario accuracy)
fprintf("Maximum RMSEs\n")
fprintf("EDMD (Trajectory): %.5f\n", max(RMSEs_EDMD_traj(not_nan_traj)))
fprintf("EDMD (Uniform): %.5f\n", max(RMSEs_EDMD_uni(not_nan_uni)))
fprintf("DE: %.5f\n", max(RMSEs_DE(not_nan_DE)))
fprintf("CCK: %.5f\n\n", max(RMSEs_CCK(not_nan_CCK)))

% number of trials where each method wins lowest RMSE
% wins(1) = number of trials where EDMD_traj got the lowest RMSE
% wins(2) = ditto for EDMD_uni
% wins(3) = ditto for DE
% wins(4) = ditto for CCK
wins = zeros(4, 1);

RMSEs = [RMSEs_EDMD_traj; RMSEs_EDMD_uni; RMSEs_DE; RMSEs_CCK];
[~, idx] = min(RMSEs, [], 1);

for i = 1:4
    wins(i) = sum(idx == i);
end

fprintf("Number of trials where each method wins with lowest RMSE\n")
fprintf("EDMD (Trajectory): %d\n", wins(1))
fprintf("EDMD (Uniform): %d\n", wins(2))
fprintf("DE: %d\n", wins(3))
fprintf("CCK: %d\n", wins(4))

if (saveRMSEs)
    save RMSEs.mat RMSEs
end

runtime = toc;
fprintf("Total runtime for mainCartPend.m: %.3f seconds\n", runtime)