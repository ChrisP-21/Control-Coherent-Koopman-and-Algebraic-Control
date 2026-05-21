% CT_EDMD.m

% uses state, input, state derivative, and output data to compute (A, B, C)
% using gEDMD

% Inputs
% X: state data
% U: input data
% X_dot: state derivative
% Y: output data
% psi: vector of lifting functions of the state variables only. Represented
% as a column of symbolic expressions
% x: symbolic column of state variables

% Outputs
% A: state matrix
% B: input matrix
% C: output matrix

function [A, B, C] = CT_EDMD(X, U, X_dot, Y, psi, x)
    disp("Computing Koopman Model using EDMD")

    N = length(psi);
    m = size(U, 1);

    % assert that all the data matrices are in appropriate dimensions
    Nsamp = size(X, 2);
    if (Nsamp ~= size(U, 2) || Nsamp ~= size(X_dot, 2) || Nsamp ~= size(Y, 2))
        error("X, U, X_dot, and Y must have the same number of columns");
    end

    if (size(X, 1) ~= size(X_dot, 1))
        error("X and X_dot must have the same number of rows");
    end

    % converting the lifting function vector and its Jacobian Dpsi into
    % function handles
    Psi = matlabFunction(psi, 'Vars', {x});

    Dpsi = jacobian(psi, x);
    Dpsi_func = matlabFunction(Dpsi, 'Vars', {x});

    % these two matrices are used for computing A and B
    psiXU_psiXU = 0;
    psidX_psiX = 0;

    % these two matrices are used for computing C
    Y_psiX = 0;
    psiX_psiX = 0;

    for i = 1:Nsamp
        % evaluate the lifting function vector and its Jacobian at each
        % state data point
        psiX = Psi(X(:, i));
        DpsiX = Dpsi_func(X(:, i)) * X_dot(:, i);

        psiXU_psiXU = psiXU_psiXU + [psiX; U(:, i)] * [psiX; U(:, i)].'/Nsamp;
        psidX_psiX = psidX_psiX + [DpsiX; zeros(m, 1)] * [psiX; U(:, i)].'/Nsamp;

        Y_psiX = Y_psiX + Y(:, i) * psiX.'/Nsamp;
        psiX_psiX = psiX_psiX + psiX * psiX.'/Nsamp;
    end

    % compute A, B, and C matrices
    K = psidX_psiX / psiXU_psiXU;

    A = K(1:N, 1:N);
    B = K(1:N, (N + 1):end);
    C = Y_psiX / psiX_psiX;
end