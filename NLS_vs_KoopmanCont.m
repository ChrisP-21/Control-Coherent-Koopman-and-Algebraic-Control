% NLS_vs_Koopman.m
% compares the Koopman model represented by KSS to the nonlinear model
% represented by NLS

% NLS = NonLinear System model
% NLS{1} = symbolic expression of f(x, u)
% NLS{2} = symbolic expression of h(x)

% KSS = Koopman State-Space model
% KSS{1} = A, KSS{2} = B, KSS{3} = C

% psi = symbolic vector of lifting functions
% x = symbolic vector of state variables
% u = symbolic vector of input variables
% dt = sampling period in seconds
% TF = final time in seconds

% x0_test = 2D array storing initial conditions for the sample trials
% (number of columns = number of trials)
% u_test = cell array of test input waveforms, each expressed as a
% function handle with time t as the input

% pltset = plot settings (cell array)
% pltset{1} = 0 or 1; 1 means plot, 0 means no plot
% pltset{2} = cell array of plot titles
function RMSE_array = NLS_vs_KoopmanCont(NLS, KSS, psi, x, u, dt, TF, x0_test, u_test, pltset)
    disp("Comparing Nonlinear Model to Koopman Model")

    % make function handles for psi, f, and h
    psi_func = matlabFunction(psi, 'Vars', {x});
    f = matlabFunction(NLS{1}, 'Vars', {x, u});
    h = matlabFunction(NLS{2}, 'Vars', {x});

    % extract A, B, C matrices from the Koopman model
    A = KSS{1};
    B = KSS{2};
    C = KSS{3};

    % computing a 2nd C matrix to obtain all the states (in order to
    % properly propagate the states in the simulations)
    bounds = [-inf(length([x; u]), 1), inf(length([x; u]), 1)];
    S = MDI_TP(x * psi', bounds, [x; u], "Gauss-Hermite", 20);
    R = MDI_TP(psi * psi', bounds, [x; u], "Gauss-Hermite", 20);
    Cx = S / R;

    t = 0:dt:TF;
    Ntrials = size(x0_test, 2);
    
    % array of RMSE values for each trial
    RMSE_array = zeros(1, Ntrials);

    % simulate trials
    warning('off', 'all')
    waiting = waitbar(0, "Starting...");
    for i = 1:Ntrials
        u_samp = u_test{i};

        % nonlinear model simulation
        [~, xOut] = ode45(@(t, x) f(x, u_samp(t)),  t, x0_test(:, i));
        yOut = h(xOut.');

        % Koopman model simulation
        [~, xOutPred] = ode45(@(t, x) Cx*(A*psi_func(x) + B*u_samp(t)), t, x0_test(:, i));
        yOutPred = C * psi_func(xOutPred.');

        % if Koopman model simulation is cut abruptly somehow, then fill
        % the rest of the yOutPred array with NaN values
        if (size(yOut, 2) > size(yOutPred, 2))
            yOutPred = [yOutPred, NaN(size(yOut, 1), size(yOut, 2) - size(yOutPred, 2))];
        end

        % compute RMSE for each trial
        RMSE_array(i) = sqrt(mean(sum((yOut - yOutPred).^2, 1)));
        waitbar(i/Ntrials, waiting, i + " out of " + Ntrials + " trials");
    end
    close(waiting);
    warning('on', 'all')
    
    % draw a plot to compute the true and predicted trajectories in the
    % very last trial
    if (pltset{1})
        for i = 1:length(h(x))
            figure
            plot(t, [yOut(i, :); yOutPred(i, :)])
            grid on
            xlabel("t [s]")
            ylabel("y_" + i)
            title(pltset{2}{i})
            legend(["True", "Predicted"], "Location", 'best');
        end
    end
end