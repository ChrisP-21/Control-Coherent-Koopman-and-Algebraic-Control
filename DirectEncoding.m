% directEncoding.m
% uses Direct Encoding to compute (A, B, C) matrices for the DE model of a
% continuous-time dynamical system represented by f and h

% Inputs
% f: symbolic continuous-time state dynamics expression. Assumes the
% system's state equation is of the form x_dot = f(x, u)
% h: symbolic expression of output. Assumes the system's output equation is
% of the form y = h(x)
% x: symbolic array of state variables (1D column vector)
% psi: symbolic array of lifting functions of x (1D column vector)
% u: symbolic array of input variables (1D column vector)
% integrationMethod: method of integration to compute inner products for
% the Gram and cross-Gram matrices
% Available integration methods:
% "QMC" = Quasi Monte Carlo (QMC)
% "MDI-TP" = Multilevel Dimension Iteration Tensor Product (MDI-TP)
% int_params: cell array of parameters for integration method

% Outputs
% A = state matrix
% B = input matrix
% C = output matrix
function [A, B, C] = DirectEncoding(f, h, x, psi, u, integrationMethod, int_params)
    disp("Computing Koopman Model using Direct Encoding")

    % symbolic expressions for time derivatives of psi
    f_sym = [psi; u];
    Df_sym = jacobian(f_sym, [x; u]) * [f; zeros(length(u), 1)];
    Df_sym = simplify(Df_sym);

    % compute Gram and cross-Gram matrices
    switch (integrationMethod)
        case "QMC"
            % Keep as backup
            bounds = int_params{1};
            w = int_params{2};
            P = int_params{3};

            Q = QuasiMonteCarlo(Df_sym * f_sym' * w, bounds, [x; u], P);
            R = QuasiMonteCarlo(f_sym * f_sym' * w, bounds, [x; u], P);
            S = QuasiMonteCarlo(h * psi' * w, bounds, [x; u], P);
        case "MDI-TP"
            % The go-to integration method for now
            bounds = int_params{1};
            w = int_params{2};
            rule = int_params{3};
            level = int_params{4};

            Q = MDI_TP(Df_sym * f_sym' * w, bounds, [x; u], rule, level);
            R = MDI_TP(f_sym * f_sym' * w, bounds, [x; u], rule, level);
            S = MDI_TP(h * psi' * w, bounds, [x; u], rule, level);
        otherwise
            error("Error: invalid integrationMethod value");
    end

    % obtain state-space matrices A, B, and C
    m = length(u);
    K = Q / R;
    A = K(1:end-m, 1:end-m);
    B = K(1:end-m, (end-m+1):end);
    C = S / R(1:end-m, 1:end-m);
end