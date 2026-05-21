# Control-Coherent-Koopman-and-Algebraic-Control
This repository contains MATLAB code that generates Koopman generator models through Control Coherent Koopman (CCK) modeling and combines them with Algebraic Control to form an accurate, robust control strategy for nonlinear systems, with theoretical accuracy and robustness guarantees.

The code was developed and tested in MATLAB 2024b.

Software requirements:
- MATLAB
	- Symbolic Math Toolbox
	- Parallel Computing Toolbox
	- Control System Toolbox
	- Robust Control Toolbox
	- Statistics and Machine Learning Toolbox (if you wish to use QMC instead of MDI-TP)
	- Progress bar (cli, gui, parfor) (link in References)
- Simulink

To help you navigate this project, the most important files are listed below:
- mainCartPend.m: the main file that handles the Koopman modeling and prediction accuracy tests for the cart pendulum example. This file computes models through Direct Encoding (DE), EDMD, and CCK and compares them with respect to prediction accuracy.
- testCartPendAlgCon.m: the file that designs an H-infinity controller for the CCK model of the actuator-augmented cart pendulum to then make an approximate inverse of the CCK model. Essential for testCartPendAlgConSim.slx.
- testCartPendAlgConSim.slx: the Simulink file that compares the combination of CCK and algebraic control to backstepping as both strategies are applied to the actuator-augmented cart pendulum. Go to Model Properties-->Callbacks to see the cart pendulum and controller parameter settings in InitFcn and simulation post-processing in StopFcn.

The rest of the files are not as important, but serve as helper functions for the following purposes:
- hermiteBasis.m: generates a symbolic array of multivariate Hermite polynomials.
- DirectEncoding.m: applies DE (or CCK depending on the plant definition) to a nonlinear system to produce (A, B, C) matrices for a Koopman LTI model.
- QuasiMonteCarlo.m: computes the required inner products for DE/CCK using Quasi Monte Carlo (QMC).
- MDI_TP.m: computes the required inner products for DE/CCK using Multilevel Dimension Iteration Tensor Product (MDI-TP).
- NLS_vs_KoopmanCont.m: compares a Koopman model to its original nonlinear system through several simulation trials and returns RMSE values for each trial.
- collectDataCT_Traj.m: generates state, input, state derivative, and output data through trajectory sampling.
- collectDataCT_Uni.m: generates state, input, state derivative, and output data through uniform random sampling.
- CT_EDMD.m: (continuous-time) EDMD that uses the data from either collectDataCT_Traj.m or collectDataCT_Uni.m to make a Koopman model.
- cartPendBackstep.m: designs a backstepping controller for the actuator-augmented cart pendulum.

References
HyunGwang Cho (2026). Progress bar (cli, gui, parfor) (https://www.mathworks.com/matlabcentral/fileexchange/121363-progress-bar-cli-gui-parfor), MATLAB Central File Exchange.
