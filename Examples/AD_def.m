function Example = AD_def()
% function Example = AnticipationDynamics_2ndMC_def()
% Define associated terms for Anticipation Dynamics

% System
sysInfo.name            = 'AD';
sysInfo.d               = 2;                                                                       % the dimension for the opinion (state) vecor
sysInfo.N               = 10;                                                                      % the number of agents
% energy_kind              = 3;                                                                       % choose a linear energy             
% [U_prime, U_dprime, T_f] = AD_get_energy_and_T_f(energy_kind);
% tau                      = 0.1;
% sysInfo.V_map           = {{@(state_i, state_j) AD_feature_map(state_i, state_j)}};
% sysInfo.V_dim           = 1;
sysInfo.phiE            = {@(r) AD_phiE(r)};
sysInfo.phiA            = {@(r) AD_phiA(r)};                                         % alignment based interaction
sysInfo.phi             = {@(r) AD_phiE(r), @(r) AD_phiA(r), @(v, d) AD_force(v, d)};    
sysInfo.phi_type        = 'EA';
sysInfo.K               = 1;                                                                       % number of classes
sysInfo.ode_order       = 2;                                                                       % order of the ODE system
sysInfo.type_info       = ones(1, sysInfo.N);                                                     % class function mapping agent index to it class index
sysInfo.agent_mass      = ones(sysInfo.N, 1);
% sysInfo.has_noise       = false;                                                                   % no stochastic noise
sysInfo.mu0             = @() AD_init_config(sysInfo.d, sysInfo.N);                              % distribution of initial conditions
sysInfo.T_f             = 20;                                                                      % the time for integration, t = T_f should be (most likely) for the system to reach steady state
sysInfo.flagxi          = 0;
sysInfo.domain          = [0, 10];
sysInfo.type = 1;
sysInfo.type_info = ones(sysInfo.N,1);

% ODE solver
solverInfo.time_span    = [0, sysInfo.T_f];                                                       % put it into the time_span vector, always starting from 0
solverInfo.rel_tol      = 1e-8;
solverInfo.abs_tol      = 1e-11;
solverInfo.option = odeset('RelTol',1e-5,'AbsTol',1e-6);


% Observations
obsInfo.L               = 5;                                                                     % observe (equi-spaced) times
obsInfo.M               = 2;                                                                     % # trajectories with random initial conditions for learning interaction kernel
obsInfo.M_rhoT          = 2000;                                                                    % # trajectories with random initial conditions to compute approximate \rho_T
obsInfo.T_0             = 0;                                                                       % Observation time interval [T_0, T]
obsInfo.T               = sysInfo.T_f/2;                                                          % Observation time interval [T_0, T]
obsInfo.time_vec        = linspace(obsInfo.T_0, obsInfo.T, obsInfo.L);                          % time instances at which discrete observation is made
obsInfo.use_derivative  = true;                                                                    % indicator of the availability of derivative data
% obsInfo.hist_num_bins   = 500;                                                                    % number of bins for estimating \rho^L_T for all three different interactions
obsInfo.obs_noise       = 0.1;
obsInfo.mu_trajnoise    = @(traj,sigma) trajUnifNoiseMultiplicative(traj, sigma);
% obsInfo.mu_dtrajnoise   = @(traj,sigma) trajUnifNoiseMultiplicative(traj, sigma);
obsInfo.rho_T_histedges    = linspace(0,sysInfo.domain(2),1000);  % a rather arbitrary set of bins to track estimators of \rho^L_T


% package the data
Example.sysInfo         = sysInfo;
Example.solverInfo      = solverInfo;
Example.obsInfo         = obsInfo;
end
