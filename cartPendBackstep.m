% cartPendBackstep.m
% 
% designs a backstepping controller for the actuator-augmented cart
% pendulum. Depends on backstepping constants a1, ..., a4, the reference
% trajectory and its derivatives, and, of course, the states in the 
% actuator-augmented cart pendulum.

% IMPORTANT NOTE: ONLY RUN THIS FILE IF THE BACKSTEPPING CONTROLLER BLOCK
% IS MISSING IN testCartPendAlgConSim.slx!!!!

%%

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

syms x1 x2 u real
x = [x1; x2];

syms p1 p2 real
p = [p1; p2];

% state equation
b = -100; % low-pass filter pole
f_cont2 = [p2;
           -b^2*p1 + 2*b*p2 + b^2*u;
           x2;
           (-k1*sin(x1) - k2*x2 - cos(x1)/k3 .* (p1 + k4*sin(x1).*x2.^2))./(1 - k5*cos(x1).^2)];

% output equation
y = x1;

%%

% designing a backstepping controller

% constants in backstepping controller
syms a [1 4] positive

% reference trajectory and its derivatives
syms r r1d r2d r3d r4d real

% backstepping "states"
syms z [1 4] real;

rnd = [r, r1d, r2d, r3d, r4d];

% z variables expressed in terms of anything but z variables
z_def = [r - y, 0, 0, 0];

% Lyapunov function candidate
V = 0.5*z(1)^2;

% define z2
dV = jacobian(V, z(1))*jacobian(z_def(1), [rnd(1:4), p.', x.'])*[rnd(2:5).'; f_cont2];
z_def(2) = a(1)*z(1) + jacobian(subs(z_def(1), z, z_def), [rnd(1:4), y])*[rnd(2:5).'; jacobian(y, [p; x])*f_cont2];
z_def(2) = subs(z_def(2), z, z_def);

% desired variable (set it to u to get the backstepping controller)
desired_var = u;

i = 2;
while (~has(dV, desired_var) && i < 5)
    % modify the Lyapunov function candidate and compute its time
    % derivative
    V = V + 0.5*z(i)^2;
    dV = subs(jacobian(V, z))*jacobian(z_def.', [rnd(1:4), p.', x.'])*[rnd(2:5).'; f_cont2];

    if (has(dV, desired_var))
        % solve for u
        temp = coeffs(dV, z(i));
        temp = temp(end);
        eqn = temp == -a(i)*z(i);
        u_con = solve(eqn, desired_var);
        u_con = subs(u_con, z, z_def);
    else
        % define the next z variable
        z_def(i + 1) = a(i)*z(i) + z(i - 1) + jacobian(z_def(i), [rnd(1:4), p.', x.'])*[rnd(2:5).'; f_cont2];
        z_def(i + 1) = subs(z_def(i + 1), z, z_def);
    end

    i = i + 1;
end

u_con = simplify(u_con);

%%

% making a MATLAB function block in Simulink that implements the
% backstepping controller
matlabFunctionBlock('testCartPendAlgConSim/CartPendBackstepping', u_con, 'Vars', {[p; x], r, r1d, r2d, r3d, r4d, a1, a2, a3, a4});