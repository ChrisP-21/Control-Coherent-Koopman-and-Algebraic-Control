clc
clear
close all

% testCartPendAlgCon.m

% designs an algebraic controller for the actuator-augmented cart pendulum

%%
if (exist("K_CCK.mat", 'file'))
    load("K_CCK.mat");

    A = A_CCK;
    B = B_CCK;
    C = C_CCK(1, :);

    % state-space model of nominal plant
    Pss = ss(A, B, C, 0);

    psi_func = matlabFunction(psi_CCK, 'Vars', {[p; x]});
    
    % (Sensitivity) Performance Weight Filter WP
    Mp = 3; % upper limit for water bad effect
    ep = 1e-3; % for getting the lowest dB -- approximate integration
    wb = 2*pi*0.5; % cut off frequency
    k = 1; % order of approximate integration (increase for steeper Bode plot)
    s = tf('s');
    WP = ss((((s/Mp^(1/k))+wb)/(s+(wb*ep^(1/k))))^k);

    % synthesize H-infinity controller
    [K,CL,GAM] = mixsyn(Pss,WP,[],[],'DISPLAY','on','METHOD','RIC');

    % just for plotting
    % If S is under 1/WP, then K stabilizes Pss
    L = ((Pss)*K); 
    S = 1/(1+L);
    T = 1 - S;
    figure;bodemag(T)
    hold on
    bodemag(S)
    bodemag(1/WP)
    legend(["T", "S", "1/W_P"]);

    loops = loopsens(Pss,K); % use this to compute sensitivity functions

    % approximate stable inverse
    Xi_a = loops.Si*K;
else
    % run mainCartPend.m if K_CCK.mat can't be found
    mainCartPend;
    save K_CCK.mat A_CCK B_CCK C_CCK f_cont2 h psi_CCK p x u
    testCartPendAlgCon;
end