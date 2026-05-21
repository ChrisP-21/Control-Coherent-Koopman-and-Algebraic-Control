% collectDataCT_Traj.m

% takes uniformly randomly selected initial conditions inside a
% hypercube whose bounds are specified in x_min, x_max, u_min, u_max and
% simulates a dynamical system with these initial conditions to produce
% more data. Data for the state, input, state derivative, and output are
% recorded afterwards.

% Inputs
% f_cont: continuous-time dynamics function in the equation x_dot = f(x,
% u). Must be a column of symbolic expressions
% h: output function in y = h(x). Also a column of symbolic expressions
% x: symbolic column vector of state variables
% u: symbolic column vector of input variables
% Ntraj: number of trajectories
% Nsamp: number of sampling points per trajectory
% dt: time step size
% x_min: column vector of minimum values for each state variable (ith row
% corresponds to ith state variable)
% x_max: column vector of maximum values for each state variable
% u_min: column vector of minimum values for each input variable (ith row
% corresponds to ith input variable)
% u_max: column vector of maximum values for each input variable

% Outputs
% X: matrix of sampling points where X(i, j) = value of ith state variable
% at jth sampling point
% U: matrix of sampling points where U(i, j) = value of ith input variable
% at jth sampling point
% X_dot: matrix of sampling points where X_dot(i, j) = value of ith state 
% variable derivative at jth sampling point
% Y: matrix of sampling points where Y(i, j) = value of ith output variable
% at jth sampling point
function [X, U, X_dot, Y] = collectDataCT_Traj(f_cont, h, x, u, Ntraj, Nsamp, dt, x_min, x_max, u_min, u_max)
    disp("Collecting trajectory data for EDMD")

    % assert that the state and input bound vectors have appropriate
    % dimensions
    n = length(x);
    m = length(u);
    if (n ~= size(x_min, 1) || n ~= size(x_max, 1))
        error("The number of rows in x_min and x_max must match the number of state variables in vector x")
    end
    if (m ~= size(u_min, 1) || m ~= size(u_max, 1))
        error("The number of rows in u_min and u_max must match the number of state variables in vector u")
    end

    % convert continuous-time dynamics function f_cont and output function
    % h into function handles
    f = matlabFunction(f_cont, 'Vars', {x, u});
    h_func = matlabFunction(h, 'Vars', {x});

    % randomly select initial conditions and randomly generate input
    % waveform for training
    t_samp = dt * (0:(Nsamp - 1));
    x0_rand = x_min + (x_max - x_min) .* rand(n, Ntraj);
    u_samp = rand(m, Nsamp);
    u_randi = cell(m, 1);
    for i = 1:m
        u_randi{i} = @(t) interp1(t_samp, u_samp(i, :), t);
    end
    u_rand = @(t) cellfun(@(func) func(t), u_randi);

    % collect state and input data via ode45 simulation
    X = zeros(n, Ntraj*Nsamp);
    U = zeros(m, Ntraj*Nsamp);
    for i = 1:Ntraj
        U(:, (1 + (i - 1)*Nsamp):(i*Nsamp)) = u_samp;

        [~, xOut] = ode45(@(t, x) f(x, u_rand(t)), t_samp, x0_rand(:, i));
        xOut = xOut.';

        if (size(xOut, 2) < Nsamp)
            xOut = [xOut, NaN(n, Nsamp - size(xOut, 2))];
        end

        X(:, (1 + (i - 1)*Nsamp):(i*Nsamp)) = xOut;
    end

    % evaluate state derivative and output at each state and input pair in
    % (X, U)
    X_dot = zeros(size(X));
    Y = zeros(length(h), Ntraj*Nsamp);
    for i = 1:(Ntraj*Nsamp)
        X_dot(:, i) = f(X(:, i), U(:, i));
        Y(:, i) = h_func(X(:, i));
    end
end