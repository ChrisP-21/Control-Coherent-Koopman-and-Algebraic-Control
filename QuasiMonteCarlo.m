% QuasiMonteCarlo.m
% element-wise integration using Quasi Monte Carlo integration.

% Inputs
% integrand: symbolic expression of the integrand. Can be scalar, vector,
% or matrix
% bounds: 2D array representing lower and upper integral bounds
% bounds(i, 1) = lower bound for ith variable
% bounds(i, 2) = upper bound for ith variable
% vars: symbolic array of variables (1D column vector)
% P: number of sampling points
% disp: boolean variable to decide whether to display progress

% Outputs
% I = integral values

function I = QuasiMonteCarlo(integrand, bounds, vars, P)
    d = length(vars);

    % generate sampling points
    if (d > 6)
        X = sobolset(d);
    else
        X = haltonset(d);
    end

    X = X(2:(P + 1), :);

    % weights
    W = ones(P, 1);

    % transform sampling points based on integration bounds
    for i = 1:d
        a = bounds(i, 1);
        b = bounds(i, 2);
        A = isinf(a);
        B = isinf(b);

        % both bounds are infinite
        if (A && B)
            W = W .* (pi * sec(pi/2 * (2*X(:, i) - 1)).^2);
            X(:, i) = tan(pi/2 * (2*X(:, i) - 1));

        % lower bound is finite, upper bound is infinite
        elseif (~A && B)
            W = W .* (1 ./ (1 - X(:, i)).^2);
            X(:, i) = a + 1 ./ (1 - X(:, i)) - 1;
            
        % lower bound is infinite, upper bound is finite
        elseif (A && ~B)
            W = W .* (-1 ./ (1 - X(:, i)).^2);
            X(:, i) = b - 1 ./ (1 - X(:, i)) + 1;
        
        % both bounds are finite
        else
            W = W .* (b - a);
            X(:, i) = a + (b - a)*X(:, i);
        end
    end

    % function handles for integrand
    N = numel(integrand);
    f = cell(N, 1);
    PB = ProgressBar(N, taskname="Making Function Handles for the Integrand");
    parfor i = 1:N
        f{i} = matlabFunction(integrand(i), 'Vars', {reshape(vars, 1, [])});
        count(PB);
    end

    % Computes I element-by-element
    I = zeros(size(integrand));

    PB = ProgressBar(N, taskname="Computing Integrals");
    parfor l = 1:N
        I(l) = mean(f{l}(X) .* W, 1);
        
        count(PB);
    end
end