% MDI_TP.m
% element-wise integration using MDI-TP.

% Inputs
% integrand: symbolic expression of the integrand. Can be scalar, vector,
% or matrix
% bounds: integral bounds (2D array with the same number of rows as the
% number of variables in vars. ith row corresponds to ith variable in vars)
% vars: symbolic array of integration variables (1D column vector)
% rule: quadrature rule
% Available options for the quadrature rule:
%   "Gauss-Legendre"
%   "Gauss-Hermite"
%       - Don't include the Gaussian weight for efficiency.
%       - Also, this rules assumes the Gaussian weight is exp(-vars.'*vars)
%   "Tanh-Sinh"
%       - Step size is 2^-6
%       - Number of sampling points = 2*level/(step size) + 1
%   Invalid string --> default to "Gauss-Legendre"
% level: quadrature precision level

% Outputs
% I = integral values

function I = MDI_TP(integrand, bounds, vars, rule, level)
    d = length(vars);
    if (size(bounds, 1) ~= d)
        disp("bounds must be a 2D array with the same number of rows" + ...
            " as the number of variables in vars")
        return
    end

    % a Boolean variable to decide whether to change the points and weights
    % according to the integral bounds
    changeYW = 1;

    % Gauss-Hermite
    if (rule == "Gauss-Hermite")
        % points and weights of Gauss-Hermite quadrature rule
        H_n = hermiteH(level, vars(1));
        y_i = double(root(H_n, vars(1)));

        H_n1 = hermiteH(level - 1, vars(1));
        H_n1 = matlabFunction(H_n1, "Vars", vars(1));
        w = @(X) (2^(level-1) * factorial(level) * sqrt(pi)/level^2) ./ H_n1(X).^2;
        w0_i = w(y_i).';
        
        changeYW = 0;

    % Tanh-Sinh
    elseif (rule == "Tanh-Sinh")
        h = 2^-6;
        N = round(level/h);

        kh = (-N:N)*h;
        y_i = tanh(pi/2 * sinh(kh)).';
        w0_i = (pi/2 * h * cosh(kh))./(cosh(pi/2 * sinh(kh)).^2);

    % default: Gauss-Legendre
    else
        % points and weights of Gauss-Legendre quadrature rule
        P_n = legendreP(level, vars(1));
        y_i = double(root(P_n, vars(1)));

        dP = diff(P_n);
        w0_i = double(2 ./ ((1 - y_i.^2) .* subs(dP, vars(1), y_i).^2)).';
    end

    % point and weight matrices
    % X_i(:, k) = evaluation points for kth variable
    % W_i(k, :) = weights for kth variable
    n = length(y_i);
    X = zeros(n, d);
    W = zeros(size(X.'));

    if (changeYW)
        for i = 1:d
            a = bounds(i, 1);
            b = bounds(i, 2);
            A = isinf(a);
            B = isinf(b);

            % both bounds are finite
            if (~A && ~B)
                X(:, i) = a + (b - a)/2 * (y_i + 1);
                W(i, :) = (b - a)/2 * w0_i;

            % lower bound is finite, upper bound is infinite 
            elseif (~A && B)
                X(:, i) = a + (1 + y_i)./(1 - y_i);
                W(i, :) = w0_i .* (2 ./ (1 - y_i).^2).';

            % lower bound is infinite, upper bound is finite
            elseif (A && ~B)
                X(:, i) = b - (1 + y_i)./(1 - y_i);
                W(i, :) = w0_i .* (-2 ./ (1 - y_i).^2).';

            % both bounds are infinite
            else
                X(:, i) = tan(pi/2 * y_i);
                W(i, :) = pi/2 * w0_i .* (sec(pi/2 * y_i).^2).';
            end
        end
    else
        for i = 1:d
            X(:, i) = y_i;
            W(i, :) = w0_i;
        end
    end

    vars = reshape(vars, 1, []);

    % scalar case (no parallel computing)
    if (isscalar(integrand))
        integrand_i = collect(integrand);
        
        % MDI loop
        PB = ProgressBar(d, taskname="Computing Integral");
        for k = d:-1:2
            integrand_i = collect(W(k, :) * subs(integrand_i, vars(k), X(:, k)), vars(1:(k - 1)));
            count(PB);
        end

        I = double(W(1, :) * subs(integrand_i, vars(1), X(:, 1)));
        count(PB);
        return;
    end
    
    I = zeros(size(integrand));
    N = numel(integrand);
    PB = ProgressBar(N, taskname="Computing Integrals");
    cellvars = num2cell(vars);

    % compute I element-by-element
    parfor l = 1:N
        
        integrand_i = collect(integrand(l));
        
        % MDI loop
        for k = d:-1:2
            integrand_i = collect(W(k, :) * subs(integrand_i, vars(k), X(:, k)), vars(1:(k - 1)));
        end

        I(l) = double(W(1, :) * subs(integrand_i, vars(1), X(:, 1)));

        count(PB);
    end
end