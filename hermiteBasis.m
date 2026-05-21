% hermiteBasis.m

% generates a symbolic array of multivariate Hermite polynomials

% Inputs
% state_vars: vector of state variables (1D symbolic vector)
%       Note: each element in state_vars can be a symbolic variable (e.g.
%       x1) or a symbolic expression (e.g. cos(x1), sin(x1))
% max_deg: maximum total degree for any function in funcs_sym
% pow_constraints: vector of max powers for each state variable (must have
% same number of elements as state_vars)

% Output
% funcs_sym: symbolic array of multivariate Hermite polynomials
function funcs_sym = hermiteBasis(state_vars, max_deg, pow_constraints)

n = length(state_vars);

if (length(pow_constraints) ~= n)
    error("Error: Power constraint array must have the same number " + ...
        "of elements as the state variable array.")
end

% generate 2-D array of powers to be used for generating multivariate
% Hermite polynomials
% HOT_POW(i, j) = power of jth variable in ith multivariate Hermite
% polynomial

HOT_POW = [];
for i = 1:max_deg
    HOT_POW = [HOT_POW; hot_pow(i, zeros(1, n), 1, pow_constraints)];
end

% use HOT_POW to generate Hermite polynomial array
N_pow = size(HOT_POW, 1);
funcs_sym = sym(zeros(N_pow, 1));

hermiteSym = hermiteH(0:max(pow_constraints), state_vars(1));

for i = 1:N_pow
    funcs_sym(i) = hermiteSym(HOT_POW(i, 1) + 1)/2;
    for l = 2:n
        funcs_sym(i) = funcs_sym(i) .* subs(hermiteSym(HOT_POW(i, l) + 1), state_vars(1), state_vars(l));
    end
end

% Computes all the combinations of nonnegative integers that sum to Q
% Q = desired sum for the integers in A
% A = current array of nonnegative integers
% k = starting index in A
% HOT_POW_MAX = array of maximal constraints on the elements of A
function HOT_POW = hot_pow(Q, A, k, HOT_POW_MAX)
    if (all(A <= HOT_POW_MAX) && sum(A) == Q)
        HOT_POW = A;
    elseif (all(A <= HOT_POW_MAX) && sum(A) < Q)
        HOT_POW = [];
        for j = k:length(A)
            B = A;
            B(j) = B(j) + 1;
            HOT_POW = [HOT_POW; hot_pow(Q, B, j, HOT_POW_MAX)];
        end
    else
        HOT_POW = [];
    end
end

end