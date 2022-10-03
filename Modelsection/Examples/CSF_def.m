function Example = CSF_def()
% function Example = OpinionDynamicsCont_def()
% Define associated terms for OPINION DYNAMICS Continuous

% (c)Sui.Tang, modified from code by M. Zhong

% System
sysInfo.name            = 'CSF';                                                   % name of the dynamics
sysInfo.d               = 2;                                                                       % dimension for the state vector (in this case, opinion vector)
sysInfo.N               = 10;                                                                      % # of agents
sysInfo.phi             = {@(r)CS_kernel(r), @(v,d) CS_force(v,d)};                                               % energy based interaction
sysInfo.phi_type        = 'A';
sysInfo.K               = 1;                                                                       % # of types
sysInfo.ode_order       = 2;                                                                       % order of the ODE system
sysInfo.type_info       = ones(1, sysInfo.N);                                                     % function mapping agent index to its type index
sysInfo.RE              = [];                                                                      % energy based reulation on interactoin beween agent i and agent i'
sysInfo.mu0             = @() CS_init_config(2,1, sysInfo.d, sysInfo.N,1);                           % distribution of initial conditions
sysInfo.T_f             = 20;         
sysInfo.flagxi          = 0;
% final time the system will reach steady state
sysInfo.domain          = [0, 30];
sysInfo.type = 1;
sysInfo.type_info = ones(sysInfo.N,1);
% ODE solver
solverInfo.time_span    = [0, sysInfo.T_f];                                                       % put it into the time_span vector, always starting from 0
solverInfo.option = odeset('RelTol',1e-5,'AbsTol',1e-6);

% Observations
obsInfo.L               = 3;            
obsInfo.M               = 3;        % # trajectories with random initial conditions for learning interaction kernel
obsInfo.time_vec = 0:5:10;
% Observations will be up to this time
obsInfo.use_derivative  = true;                                                                   % indicator of the availability of derivative data
obsInfo.obs_noise       = 0;
obsInfo.mu_trajnoise    = @(traj,sigma) trajUnifNoiseAdditive( traj, sigma );
obsInfo.rho_T_histedges    = linspace(0,sysInfo.domain(2),1000);  % a rather arbitrary set of bins to track estimators of \rho^L_T


% package the data
Example.sysInfo         = sysInfo;
Example.solverInfo      = solverInfo;
Example.obsInfo         = obsInfo;


end