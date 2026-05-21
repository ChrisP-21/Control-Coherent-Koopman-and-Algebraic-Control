% collectDataCT_Uni.m

% takes uniformly randomly selected sampling points inside a hypercube
% whose bounds are specified in x_min, x_max, u_min, u_max. This file then
% records the state, input, state derivative, and output at these sampling
% points.

% Inputs
% f_cont: continuous-time dynamics function in the equation x_dot = f(x,
% u). Must be a column of symbolic expressions
% h: output function in y = h(x). Also a column of symbolic expressions
% x: symbolic column vector of state variables
% u: symbolic column vector of input variables
% Nsamp: number of sampling points
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
function [X, U, X_dot, Y] = collectDataCT_Uni(f_cont, h, x, u, Nsamp, x_min, x_max, u_min, u_max)
    disp("Collecting uniformly random data for EDMD")

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
    % h into functions
    f = matlabFunction(f_cont, 'Vars', {x, u});
    h_func = matlabFunction(h, 'Vars', {x});

    % randomly select state and input vector values
    X = x_min + (x_max - x_min).*rand(n, Nsamp);
    U = u_min + (u_max - u_min).*rand(m, Nsamp);

    % evaluate state derivative and output at each state and input pair in
    % (X, U)
    X_dot = zeros(size(X));
    Y = zeros(length(h), Nsamp);
    for i = 1:Nsamp
        X_dot(:, i) = f(X(:, i), U(:, i));
        Y(:, i) = h_func(X(:, i));
    end
end