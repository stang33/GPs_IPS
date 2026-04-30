% Authored by Jinchao Feng and Sui Tang
% RunExamples for the main cases
clear all;
close all;

addpaths;

%% Set 2parameters
% if ispc, SAVE_DIR = [getenv('USERPROFILE'), '\DataAnalyses\LearningDynamics']; else, 
%SAVE_DIR = [getenv('HOME'), '/DataAnalyses/LearningDynamics']; 
%end % Please keep this fixed, simply create a symlink ~/DataAnalyses pointing wherever you like                           
VERBOSE                         = 1;                                                                % indicator to print certain output
time_stamp                      = datestr(now, 30);
if ~exist('Params','var'), Params = [];     end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end

%% Load example definitions
% System
sysInfo.name            = 'ODNS';                                                   % name of the dynamics
sysInfo.d               = 1;                                                                       % dimension for the state vector (in this case, opinion vector)
sysInfo.N               = 100;                                                                      % # of agents
sysInfo.phi             = {@(r)ODNS_influence(r, 1)};                                               % energy based interaction
sysInfo.Noption         = 1;
sysInfo.K               = 1;                                                                       % # of types
sysInfo.ode_order       = 1;                                                                       % order of the ODE system
sysInfo.type_info       = ones(1, sysInfo.N);                                                     % function mapping agent index to its type index
sysInfo.RE              = [];                                                                      % energy based reulation on interactoin beween agent i and agent i'
sysInfo.flagxi          = 0;

sysInfo.mu0             = @()ODNS_init_config(sysInfo.d, sysInfo.N, sysInfo.d);                           % distribution of initial conditions
sysInfo.T_f             = 10;                                                                      % final time the system will reach steady state
sysInfo.domain          = [0, 10];
sysInfo.type = 1;
sysInfo.type_info = ones(sysInfo.N,1);
% ODE solver
solverInfo.time_span    = [0, sysInfo.T_f];                                                       % put it into the time_span vector, always starting from 0
solverInfo.option = odeset('RelTol',1e-5,'AbsTol',1e-6);

% Observations
obsInfo.M               = 1;                         % chosen from 1,...,10                                              % # trajectories with random initial conditions for learning interaction kernel
obsInfo.time_vec = 0:1:5;

% Observations will be up to this time
obsInfo.use_derivative  = true;                                                                   % indicator of the availability of derivative data
% obsInfo.obs_noise       = 0;
obsInfo.mu_trajnoise    = @(traj,sigma) trajUnifNoiseAdditive( traj, sigma );
obsInfo.rho_T_histedges    = linspace(0,sysInfo.domain(2),1000);  % a rather arbitrary set of bins to track estimators of \rho^L_T
                                            % move n to learn_info
obsInfo.VERBOSE                = VERBOSE;
obsInfo.SAVE_DIR               = SAVE_DIR;

%% Learning setup (named baseLearnInfo so parfor can treat 'learnInfo' as a temporary)
baseLearnInfo.N = sysInfo.N;
baseLearnInfo.d = sysInfo.d;
baseLearnInfo.order = sysInfo.ode_order;
baseLearnInfo.name = sysInfo.name;
baseLearnInfo.jitter = 1e-6;
baseLearnInfo.dtraj ='false'; % do not compute dtraj during trajectory computation

saveON= 0;
plotON = 0;

% %% generate rhoT
% obsInfo.MrhoT = 2000;
% learnInfo.rhoLT = Generate_rhoT(sysInfo,obsInfo,solverInfo,saveON,plotON);% empirical pairwise distance

%% construct basis function
nbasis = 100;
rmax = sysInfo.domain(2)-sysInfo.domain(1);

rng(1);

% basis information
basis_info.type = 'piecewise';

switch basis_info.type
    case 'piecewise'
        edges = linspace(0,rmax,nbasis+1);
        basis_funs = localBasisFn(edges,0);
     
    % case 'harmonic'
    %     basis_funs = cell(1,nbasis);
    % 
    %      for k = 1:nbasis
    %         basis_funs{k} = @(x) cos((k-1).*x).^2;
    %      end
    % 
    % case 'poly'
    %     basis_funs = cell(1,nbasis);
    %     for k = 1:nbasis
    %         basis_funs{k} = @(x) x.^(k-1);
    %     end
    % 
    % case 'invpoly'
    %     basis_funs = cell(1,nbasis);
    %     for k = 1:nbasis
    %         basis_funs{k} = @(x) x.^(-k);
    %     end
    % 
    % case 'exp'
    %     basis_funs = cell(1,nbasis);
    %     for k = 1:ceil(nbasis)
    %         basis_funs{k} = @(x) exp(-k.*x);
    %     end
    %     % for k = 1:floor(nbasis/2)
    %     %     basis_funs{ceil(nbasis/2)+k} = @(x) exp(-x./(k+1));
    %     % end
    % case 'all'
    %     edges = linspace(0,rmax,nbasis+1);
    %     basis_funs1 = localBasisFn(edges,0);  
    % 
    %     basis_funs2 = cell(1,nbasis);
    %     for k = 1:nbasis
    %         basis_funs2{k} = @(x) x.^(k);
    %     end
    % 
    %     basis_funs3 = cell(1,nbasis);
    %     for k = 1:nbasis
    %         basis_funs3{k} = @(x) x.^(-k);
    %     end
    % 
    %     basis_funs4 = cell(1,nbasis);
    %     for k = 1:ceil(nbasis)
    %         basis_funs4{k} = @(x) exp(-k.*x);
    %     end
    % 
    %     basis_funs = [basis_funs1 basis_funs2 basis_funs3 basis_funs4];
    % 
    %     nbasis = size(basis_funs,2);
    % 
    % case 'CSkernel'
    %      basis_funs = cell(1,nbasis);
    % 
    %      for k = 1:nbasis
    %         basis_funs{k} = @(x) 1./(1 + x.^2).^(1/k);
    %      end
    % 
    % case 'Chebyshev'
    %     basis_funs = cell(1,nbasis);
    %     chebyshevpoly{1} = @(x) 1;
    %     chebyshevpoly{2} = @(x) x;
    %     chebyshevpoly{3} = @(x) 2*x.^2-1;
    %     chebyshevpoly{4} = @(x) 4*x.^3-3*x;
    %     chebyshevpoly{5} = @(x) 8*x.^4-8*x.^2+1;
    %     chebyshevpoly{6} = @(x) 16*x.^5-20*x.^3+5*x;
    %     chebyshevpoly{7} = @(x) 32*x.^6-48*x.^4+18*x.^2-1;
    %     chebyshevpoly{8} = @(x) 64*x.^7-112*x.^5+56*x.^3-7*x;
    %     chebyshevpoly{9} = @(x) 128*x.^8-256*x.^6+160*x.^4-32*x.^2+1;
    %     chebyshevpoly{10} = @(x) 256*x.^9-576*x.^7+432*x.^5-120*x.^3+9*x;
    %     chebyshevpoly{11} = @(x) 512*x.^10-1280*x.^8+1120*x.^6-400*x.^4+50*x.^2-1;
    % 
    %     if nbasis>10
    %         for k = 1:nbasis
    %             basis_funs{k} = @(x) chebyshevT(k,x);
    %         end
    %     else
    %         for k = 1:nbasis
    %             basis_funs{k} = @(x) chebyshevpoly{k}(x);
    %         end
    %     end
         
    % case 'legendre'
    %     basis_funs = cell(1,nbasis);
    %     for k = 1:nbasis
    %         basis_funs{k} = @(x) legendreP(k-1,x/rmax*2-1);
    %     end

    % case 'random_feature'
    %     basis_info.feature_type = 'gaussian';
    %     basis_info.normalize = true;
    %     basis_info.sigma = 0.1;
    % 
    %     basis = construct_random_feature_basis(rmax, nbasis, basis_info);
    %     basis_funs = basis.f;
    % 
    % case 'Gaussian'
    %     basis_funs = cell(1,nbasis);
    %     edges = linspace(0,rmax,nbasis+1);
    %     sigma = 0.1;
    %     for ind = 1 : nbasis
    %         basis_funs{ind} = @(r) exp(-(r - (edges(ind)+edges(ind+1))/2).^2 ./ (2*sigma^2));
    %     end
    % 
    % case 'nonparametric'
    %     basis_info.kernel_type = 'gaussian';
    %     r_X = [];
    %     for m = 1:learnInfo.M
    %         for l = 1:learnInfo.L
    %             for i = 1:learnInfo.N
    %                r_X = [r_X; xpath_train([1:sysInfo.d*(i-1), sysInfo.d*i+1:sysInfo.d*learnInfo.N],l,m) - repmat(xpath_train((i-1)*learnInfo.d+1:i*learnInfo.d,l,m),learnInfo.N-1,1)]; 
    %             end
    %         end
    %     end
    %     r_X = reshape(r_X,sysInfo.d,[]);
    %     r_X_norm = sqrt(sum(r_X.^2,1));
    %     basis_info.sigma = 1;
    %     basis = construct_nonparametric_basis(unique(r_X_norm), nbasis, basis_info);
    %     basis_funs = basis.f;
        
end

%% --- M sweep + parfor over trials ---
for M_curr = 1:10
    obsInfo.M = M_curr;
    fprintf('\n========== ODNS1 M = %d ==========\n', M_curr);

for nsr   = [0 0.25 0.5 0.75 1]

    %%
    nT = 10;     %number of trials
    % errorphis = zeros(2,nT);      %store errors of phis in L-infinity and L2rhoT norms
    errortrajs_train = zeros(4,nT);     %store mean and std of trajectory error in training data
    errortrajs_test = zeros(4,nT);      %store mean and std of trajectory error in testing data

    F_trials = cell(1,nT);

    obsInfo.obs_noise_nsr = nsr;   % hoisted out of parfor: needed by save() below

    parfor n = 1:nT
        learnInfo = baseLearnInfo;   % per-iteration local copy of the broadcast setup
        %% generating data
        [dxpath_test,xpath_test,dxpath_train, xpath_train]=Generate_training_data(sysInfo,obsInfo,solverInfo);

        learnInfo.rho_emp = rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);% empirical pairwise distance

        % noise = 0;
        % obsInfo.obs_noise = noise;
        % noise_to_signal_ratio = 0;
    
        dxpath_true = dxpath_train;
        % xpath_train= trajUnifNoiseAdditive(xpath_train, obsInfo.obs_noise);
        % xpath_train= trajUnifNoiseMultiplicative(xpath_train, obsInfo.obs_noise);
    
        if sysInfo.ode_order == 1
            dxpath_train = trajUnifNoiseAdditive(dxpath_train, mean(abs(dxpath_train),'all')*obsInfo.obs_noise_nsr);
        elseif sysInfo.ode_order == 2
            dxpath_train = trajUnifNoiseAdditive(dxpath_train, mean(abs(dxpath_train(sysInfo.d*sysInfo.N+1:end,:,:)),'all')*obsInfo.obs_noise_nsr);
        end
    
        fprintf('\n Done! begin to learn ......\n');
        
        X_train=[];
        Y_train=[];
        
        %  For second order system:     X_train consists of [x_1,\cdots,x_N,v_1,\cdots,v_N]
        %                               Y_train consistss of [\dot v_1,..., \dot v_N]
        
        for i = 1:obsInfo.M
            for j = 1:size(xpath_train,2)
                X_train = [X_train;xpath_train(:,j,i)];
                if sysInfo.ode_order ==2
                    Y_train = [Y_train;dxpath_train(sysInfo.d*sysInfo.N+1:end,j,i)];
                else
                    Y_train = [Y_train;dxpath_train(:,j,i)];
                    Y_true = [Y_train;dxpath_true(:,j,i)];
                end
            end
        end
        
        learnInfo.L= size(xpath_train,2);
        learnInfo.M = size(xpath_train,3);
        learnInfo.X = X_train;
        learnInfo.Y = Y_train;
        
        learnInfo.xpath_train = xpath_train;
        learnInfo.dxpath_train = dxpath_train;
        
        
        % %% visualize the training trajectory
        % 
        % if sysInfo.d ==1
        %     visualize_training_trajs_1D(sysInfo,obsInfo,solverInfo,learnInfo);
        % 
        % elseif sysInfo.d ==2
        %         visualize_training_trajs_2D(sysInfo,obsInfo,solverInfo,learnInfo);
        % 
        % end
    
    
        %% construct F
        % F = zeros(nbasis,learnInfo.d*learnInfo.N*learnInfo.L*learnInfo.M);
        F_all = [];
        F_true = [];
        % B_all = [];
        r_rho = [];
        
        for m = 1:learnInfo.M
            for l = 1:learnInfo.L
                for i = 1:learnInfo.N
                    if sysInfo.Noption == 1  %normalization including phi_i(0)
                        r_temp = xpath_train(1:sysInfo.d*learnInfo.N,l,m) - repmat(xpath_train((i-1)*learnInfo.d+1:i*learnInfo.d,l,m),learnInfo.N,1); 
                    else
                        r_temp = xpath_train([1:sysInfo.d*(i-1), sysInfo.d*i+1:sysInfo.d*learnInfo.N],l,m) - repmat(xpath_train((i-1)*learnInfo.d+1:i*learnInfo.d,l,m),learnInfo.N-1,1); 
                    end
                   r_temp = reshape(r_temp,sysInfo.d,[]);
                   r_temp_norm = sqrt(sum(r_temp.^2,1));
                   r_rho = [r_rho r_temp_norm];
                   if sysInfo.ode_order ==2
                       if sysInfo.Noption == 1
                           v_temp = xpath_train(sysInfo.d*sysInfo.N+1:sysInfo.d*sysInfo.N+sysInfo.d*learnInfo.N,l,m) - repmat(xpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m),learnInfo.N,1); 
                       else
                           v_temp = xpath_train([sysInfo.d*sysInfo.N+1:sysInfo.d*sysInfo.N+sysInfo.d*(i-1),sysInfo.d*sysInfo.N+sysInfo.d*i+1:sysInfo.d*sysInfo.N+sysInfo.d*learnInfo.N],l,m) - repmat(xpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m),learnInfo.N-1,1); 
                       end
                       v_temp = reshape(v_temp,sysInfo.d,[]);
                   end
                   F_col = zeros(nbasis,sysInfo.d);
                   F_true_col = zeros(nbasis,sysInfo.d);
                   % B_col = zeros(nbasis,1);
        
                   if strcmp('ODNS',learnInfo.name)
                       for k = 1:nbasis
                           F_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_train((i-1)*learnInfo.d+1:i*learnInfo.d,l,m)-r_temp,sysInfo.d,[]),2)';
                           F_true_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_true((i-1)*learnInfo.d+1:i*learnInfo.d,l,m)-r_temp,sysInfo.d,[]),2)';
                           % B_col(k,:) = sum(basis_funs{k}(r_temp_norm));
                       end
                   elseif strcmp('LJ',learnInfo.name)
                       for k = 1:nbasis
                           F_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_train((i-1)*learnInfo.d+1:i*learnInfo.d,l,m)-r_temp,sysInfo.d,[]),2)';
                           F_true_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_true((i-1)*learnInfo.d+1:i*learnInfo.d,l,m)-r_temp,sysInfo.d,[]),2)';
                           % B_col(k,:) = sum(basis_funs{k}(r_temp_norm));
                       end
                   elseif strcmp('CSNS',learnInfo.name)
                       for k = 1:nbasis
                           F_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m)-v_temp,sysInfo.d,[]),2)';
                           F_true_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_true(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m)-v_temp,sysInfo.d,[]),2)';
                       end
                   elseif strcmp('FMNS',learnInfo.name)
                       for k = 1:nbasis
                           F_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m)-sysInfo.phi{2}(xpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m),sysInfo.d)-r_temp,sysInfo.d,[]),2)';
                           F_true_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_true(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m)-sysInfo.phi{2}(xpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m),sysInfo.d)-r_temp,sysInfo.d,[]),2)';
                       end                   
                   end
        
                   F_all = [F_all F_col];
                   F_true = [F_true F_true_col];
                   % B_all = [B_all B_col];
                end
            end
        end
    
        F = F_all';
        F_true = F_true';
        
        normF = ones(1,nbasis);
    
        % normalizing the basis/each column of the matrix F
        % for i = 1:nbasis   
        %     normF(i) = norm(F(:,i));
        %     F(:,i) = F(:,i)./normF(i);
        %     F_true(:,i) = F_true(:,i)./normF(i);
        % end
        % F(isnan(F))=0;
    
        F_trials{n} = F_all';
    
    %% ============================================================
    % implicit-SINDy comparison baseline (faithful homogeneous version)
    % Solve F c = 0 by finding the sparsest nonzero vector in an
    % approximate null space using alternating directions method (ADM)
    % ============================================================
    
    % run_iSINDy_compare = 0;
    % 
    % if run_iSINDy_compare
    %     iS_opts = struct();
    %     iS_opts.null_dims = 1:nbasis-1;     % candidate null-space dimensions
    %     iS_opts.lambda_grid = logspace(-4, 0, 20);  % sparsity thresholds for ADM
    %     iS_opts.maxit = 500;
    %     iS_opts.tol = 1e-6;
    %     iS_opts.verbose = false;
    % 
    %     % same train/test split idea as your Bayesian model-selection stage
    %     rn = randperm(size(F,1));
    %     train_id = rn(1:round(3/4*size(F,1)));
    %     test_id  = rn(round(3/4*size(F,1))+1:end);
    % 
    %     F_train = F(train_id, :);
    %     F_test  = F(test_id,  :);
    % 
    %     [coef_iS, iS_info] = run_implicit_sindy_homogeneous(F_train, F_test, iS_opts);
    % 
    %     % convert coefficient vector into kernel estimator
    %     phi_est_iS = @(r) 0;
    %     for k = 1:nbasis
    %         if abs(coef_iS(k)) > 0
    %             phi_est_iS = @(r) phi_est_iS(r) + coef_iS(k) ./ normF(k) .* basis_funs{k}(r);
    %         end
    %     end
    % 
    %     fprintf('\nimplicit-SINDy comparison finished.\n');
    %     fprintf('  chosen null dimension = %d\n', iS_info.best_null_dim);
    %     fprintf('  chosen lambda         = %.4e\n', iS_info.best_lambda);
    %     fprintf('  validation residual   = %.4e\n', iS_info.best_residual);
    %     fprintf('  nnz(coefficients)     = %d\n', nnz(abs(coef_iS) > 1e-10));
    % end
    
        %% Sparse Bayesian regression
        % set i-th coefficient to be one
        method1 = 'FastLaplace'; 
        MScriterion1 = ones(1,nbasis)*10^10;
    
        for i = 1:nbasis
            eta = F(:, i);
            sigma = 0;
            eta = eta + sigma*randn(size(eta,1),size(eta,2));
            % Phi = -F(:,[1:i-1,i+1:nbasis]);
            rn = randperm(size(F,1));
            Phi_train = -F(rn(1:round(3/4*size(F,1))),[1:i-1,i+1:nbasis]);
            eta_train = eta(rn(1:round(3/4*size(F,1))));
            Phi_test = -F(rn(round(3/4*size(F,1))):end,[1:i-1,i+1:nbasis]);
            eta_test = eta(rn(round(3/4*size(F,1))):end);
            if sum(eta_train~=0) < 0.5*size(eta_train,1)
                continue
            end
    
            % ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
            ind = [1:i-1,i+1:nbasis];
    
            try 
                sigma2 = var(eta)/1e3;
                delta_La = 1e-10;
                lambda_init = [];
                sparsity = [];
                [weights,used_ids,loglik,sigma2,errbars,basis,selected,deleted,alpha,lambda] = FastLaplace(Phi_train,eta_train, sigma2, delta_La,lambda_init,sparsity);
                % deleted_ids = (abs(weights)./norm(weights)<threshold);
                % used_ids(deleted_ids) = [];
                % weights(deleted_ids) = [];
                % alpha(deleted_ids) = [];
    
                threshold = [];
                %threshold = 1e-3;
                %threshold = sqrt(sigma2);
                if ~isempty(threshold)
                    %if any(abs(weights)./norm(weights)<threshold)
                    if any(abs(weights)'.*vecnorm(Phi_train(:,used_ids))/sqrt(size(Phi,1))<threshold)
                        deleted_ids = (abs(weights)./norm(weights)<threshold);
                        used_ids(deleted_ids) = [];
                        [weights,used_ids_t,loglik,sigma2,errbars,basis,selected,deleted,alpha,lambda] = FastLaplace(Phi_train(:,used_ids),eta_train, sigma2, delta_La,lambda_init,sparsity);
                        used_ids = used_ids(used_ids_t);
                    end
                end
    
                temp = zeros(nbasis-1,1);
                temp(used_ids) = weights;
                xhat = temp;
    
                U = chol(Phi_train(:,used_ids)'*Phi_train(:,used_ids)./sigma2 + diag(alpha));
                Ui = inv(U);
                SIGMA = Ui * Ui';
                errbars = zeros(length(ind),1);
                errbars(used_ids) = sqrt(diag(SIGMA));
            catch
                continue
            end
    
            if sum(Phi_train(:,used_ids)~=0,"all") < nbasis
                continue
            end
    
            MScriterion1(1,i) = (sum(errbars(used_ids).^2)+sigma2)./(sum(weights.^2)+1);
        end
        
    
        [~,L_i] = min(MScriterion1(1,:)); % find the best 'i-th' basis to set its coefficient to be one
    
        %% SBR with the coefficient of i-th basis set to be one
        iA = L_i;
        eta = F(:, iA);
        sigma = 0; % add additional model noise
        eta = eta + sigma*randn(size(eta,1),size(eta,2));
        Phi = -F(:,[1:iA-1,iA+1:nbasis]);
        ibasis_funs = basis_funs(1,[1:iA-1,iA+1:nbasis]);
        ind = [1:iA-1,iA+1:nbasis];
    
        sigma2 = var(eta)/1e3 ;
        delta_La = 1e-10;
        lambda_init = [];
        sparsity = [];
        tic
        [weights,used_ids,loglik,sigma2,errbars,basis,selected,deleted,alpha,lambda] = FastLaplace(Phi,eta, sigma2, delta_La,lambda_init,sparsity);
        time_SBL = toc;
    
        threshold = [];
        if ~isempty(threshold)
            if any(abs(weights)'.*vecnorm(Phi_train(:,used_ids))/sqrt(size(Phi,1))<threshold)
                deleted_ids = (abs(weights)./norm(weights)<threshold);
                used_ids(deleted_ids) = [];
                [weights,used_ids_t,loglik,sigma2,errbars,basis,selected,deleted,alpha,lambda] = FastLaplace(Phi(:,used_ids),eta, sigma2, delta_La,lambda_init,sparsity);
                used_ids = used_ids(used_ids_t);
            end
        end
    
        temp = zeros(nbasis-1,1);
        temp(used_ids) = weights;
        xhat = temp;
    
        U = chol(Phi(:,used_ids)'*Phi(:,used_ids)./sigma2 + diag(alpha));
        Ui = inv(U);
        SIGMA = Ui * Ui';
        errbars_L = zeros(length(ind),1);
        errbars_L(used_ids) = sqrt(diag(SIGMA));
    
        est_L = xhat;
        phi_est_L = @(r) 1/normF(iA).*basis_funs{iA}(r);
        for k = 1:nbasis-1
            if est_L(k) ~= 0
                phi_est_L = @(r) phi_est_L(r) + est_L(k)./normF(ind(k)).*ibasis_funs{k}(r);
            end
        end

        coef_L = zeros(nbasis,1);
        coef_L(1:iA-1,1) = est_L(1:iA-1)./normF(1:iA-1)';
        coef_L(iA,1) = 1 / normF(iA);
        coef_L(iA+1:nbasis,1) = est_L(iA:nbasis-1)./normF(iA+1:nbasis)';
    
        dr = edges(2) - edges(1);
        phi_fast = @(r) phi_piecewise_const(r, coef_L, rmax);
    
        % [errorphis(1,n),errorphis(2,n)] = errornorms_phis_NS(sysInfo,obsInfo,learnInfo,range,phi_est);
        result_train = construct_and_compute_traj_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_fast,learnInfo.xpath_train(:,1,:));
        errortrajs_train(:,n) = [result_train.train_traj_error result_train.prediction_traj_error]';
        result_test = construct_and_compute_traj_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_fast,sysInfo.mu0());
        errortrajs_test(:,n) = [result_test.train_traj_error result_test.prediction_traj_error]';
    
        % %% Bayesian Linear Model
        % method2 = 'Gaussian'; 
        % MScriterion2 = ones(1,nbasis)*10^10;
        % 
        % for i = 1:nbasis
        %     eta = F(:, i);
        %     sigma = 0;
        %     eta = eta + sigma*randn(size(eta,1),size(eta,2));
        %     rn = randperm(size(F,1));
        %     Phi_train = -F(rn(1:round(3/4*size(F,1))),[1:i-1,i+1:nbasis]);
        %     eta_train = eta(rn(1:round(3/4*size(F,1))));
        %     Phi_test = -F(rn(round(3/4*size(F,1))):end,[1:i-1,i+1:nbasis]);
        %     eta_test = eta(rn(round(3/4*size(F,1))):end);
        % 
        %     if sum(eta_train~=0) < 0.5*size(eta_train,1)
        %         continue
        %     end
        % 
        %     % ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
        %     ind = [1:i-1,i+1:nbasis];
        % 
        %     try 
        %         [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi_train, eta_train);
        % 
        %         threshold = [];
        %         if ~isempty(threshold)
        %             if any(abs(PARAMETER.Value)'.*vecnorm(Phi_train(:,PARAMETER.Relevant))/sqrt(size(Phi,1))<threshold)
        %                 deleted_ids = (abs(PARAMETER.Value)./norm(PARAMETER.Value)<threshold);
        %                 PARAMETER.Relevant(deleted_ids) = [];
        %                 [PARAMETER_t, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi(:,PARAMETER.Relevant), eta);
        %                 PARAMETER.Relevant = PARAMETER.Relevant(PARAMETER_t.Relevant);
        %                 PARAMETER.Value = PARAMETER.Value((PARAMETER_t.Relevant));
        %             end
        %         end
        % 
        % 
        %         xhat = zeros(length(ind),1);
        %         xhat(PARAMETER.Relevant) = PARAMETER.Value;
        %         U = chol(Phi_train(:,PARAMETER.Relevant)'*Phi_train(:,PARAMETER.Relevant)*HYPERPARAMETER.beta + diag(HYPERPARAMETER.Alpha));
        %         Ui = inv(U);
        %         SIGMA = Ui * Ui';
        %         errbars = zeros(length(ind),1);
        %         errbars(PARAMETER.Relevant) = sqrt(diag(SIGMA));
        %     catch
        %         continue
        %     end
        % 
        %     if sum(Phi_train(:,PARAMETER.Relevant)~=0,"all") < nbasis
        %         continue
        %     end
        % 
        %     MScriterion2(1,i) = (sum(errbars(PARAMETER.Relevant).^2)+1./HYPERPARAMETER.beta)./(sum(PARAMETER.Value.^2)+1);
        % end
        % 
        % 
        % [~,G_i] = min(MScriterion2(1,:)); % find the best 'i-th' basis to set its coefficient to be one
        % 
        % %% BLM with the coefficient of i-th basis set to be one
        % i = G_i;
        % eta = F(:, i);
        % sigma = 0; % add additional model noise
        % eta = eta + sigma*randn(size(eta,1),size(eta,2));
        % Phi = -F(:,[1:i-1,i+1:nbasis]);
        % ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
        % ind = [1:i-1,i+1:nbasis];
        % 
        % [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi, eta);
        % 
        % threshold = [];
        % if ~isempty(threshold)
        %     if any(abs(PARAMETER.Value)'.*vecnorm(Phi_train(:,PARAMETER.Relevant))/sqrt(size(Phi,1))<threshold)
        %         deleted_ids = (abs(PARAMETER.Value)./norm(PARAMETER.Value)<threshold);
        %         PARAMETER.Relevant(deleted_ids) = [];
        %         [PARAMETER_t, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi(:,PARAMETER.Relevant), eta);
        %         PARAMETER.Relevant = PARAMETER.Relevant(PARAMETER_t.Relevant);
        %         PARAMETER.Value = PARAMETER.Value((PARAMETER_t.Relevant));
        %     end
        % end
        % 
        % xhat = zeros(length(ind),1);
        % xhat(PARAMETER.Relevant) = PARAMETER.Value;
        % 
        % U = chol(Phi(:,PARAMETER.Relevant)'*Phi(:,PARAMETER.Relevant)*HYPERPARAMETER.beta + diag(HYPERPARAMETER.Alpha));
        % Ui = inv(U);
        % SIGMA = Ui * Ui';
        % errbars_G = zeros(length(ind),1);
        % errbars_G(PARAMETER.Relevant) = sqrt(diag(SIGMA));
        % 
        % est_G = xhat;
        % phi_est_G = @(r) 1/normF(i).*basis_funs{i}(r);
        % for k = 1:nbasis-1
        %     if est_G(k) ~= 0
        %         phi_est_G = @(r) phi_est_G(r) + est_G(k)./normF(ind(k)).*ibasis_funs{k}(r);
        %     end
        % end
        % 
        % coef_G = zeros(nbasis,1);
        % coef_G(1:i-1,1) = est_G(1:i-1)./normF(1:i-1)';
        % coef_G(i,1) = 1 / normF(i);;
        % coef_G(i+1:nbasis,1) = est_G(i:nbasis-1)./normF(i+1:nbasis)';
    
    end
    
    %% save data
    out_dir = '../trajectory_data';
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    save(strcat(out_dir,'/',sysInfo.name,'1N',num2str(sysInfo.N),'d',num2str(sysInfo.d),'M',num2str(obsInfo.M),'L',num2str(length(obsInfo.time_vec)),'sigma',num2str(obsInfo.obs_noise_nsr), '.mat'),'errortrajs_train','errortrajs_test');
end
end  % end M_curr sweep
%% plot comparison of learned kernel
if run_iSINDy_compare
    % Define the range for r
    r = linspace(0, rmax, 500); % Smooth interpolation
    
    % Compute original and estimated phi values
    phi_sys = sysInfo.phi{1}(r);
    phi_L = phi_est_L(r);
    phi_G = phi_est_G(r);
    phi_iS = phi_est_iS(r);
    
    % Normalize by relative l1 norm (mean absolute value)
    % normalize_sys = 1/mean(abs(phi_sys));
    normalize_L = mean(abs(phi_sys(1:200)))/mean(abs(phi_L(1:200)));
    normalize_G = mean(abs(phi_sys(1:200)))/mean(abs(phi_G(1:200)));
    normalize_iS = mean(abs(phi_sys(1:500)))/mean(abs(phi_iS(1:500)));
    
    % phi_sys = phi_sys * normalize_sys;
    phi_L = phi_L * normalize_L;
    phi_G = phi_G * normalize_G;
    phi_iS = phi_iS * normalize_iS;
    scrsz = [1, 1, 1920, 1080];
    
    % Create figure with size matching trajectory plots
    phi_fig = figure('Name', 'Comparison of learning methods', ...
                     'NumberTitle', 'off', ...
                     'Position', [scrsz(3)/8, scrsz(4)/8, scrsz(3)*3/4, scrsz(4)*3/4]);
    
    % Set background to white
    set(gcf, 'Color', 'w');
    
    % Define improved color palette
    colors = {[0, 0, 0], [0, 0.447, 0.741], [0.85, 0.325, 0.098]};
    
    % Plot curves with enhanced styles
    p_true = plot(r, phi_sys, '-', 'LineWidth', 3, 'Color', colors{1}); % True \phi
    hold on;
    p_L = plot(r, phi_L, '--', 'LineWidth', 3, 'Color', colors{2});  % Estimated \phi (L)
    p_G = plot(r, phi_G, '-.', 'LineWidth', 3, 'Color', colors{3});  % Estimated \phi (B)
    p_iS = plot(r, phi_iS, ':', 'LineWidth', 3, 'Color', 'magenta');  % Estimated \phi (B)

    % hold off;

    % plot the density of rho
    yyaxis right                                                                                % display \rho^L_T and its estimator
    axesHandle1         = gca();
    hist1       = histogram(r_rho,100,'EdgeColor',[1 1 1],'Normalization','pdf');
    hist1.FaceAlpha = 0.3;
    
    hold off;
    
    % Apply font sizes consistent with your trajectory visualization
    set(gca, 'FontSize', 39, 'FontName', 'Helvetica', 'TickLabelInterpreter', 'latex');
    
    % Labeling with LaTeX interpreter
    yyaxis left
    xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 50, 'FontName', 'Helvetica');
    ylabel('$\phi(r)$', 'Interpreter', 'latex', 'FontSize', 50, 'FontName', 'Helvetica');
    title(['Estimation of $\phi$ (noise level = ', num2str(noise_to_signal_ratio*100), '\%)'], ...
        'Interpreter', 'latex', 'FontSize', 50, 'FontName', 'Helvetica');
    
    % Add a refined legend with SIAM-compatible font size
    legend([p_true, p_L, p_G, p_iS, hist1],{'True $\phi$', 'FastLaplace', 'Gaussian', 'implicit-SINDy', '$\rho_r$'}, ...
        'Interpreter', 'latex', 'FontSize', 33, 'Location', 'northeast');
    
    % Adjust axes and grid
    xlim([0, rmax]);
    yyaxis left
    ylim([min([phi_sys,phi_est_L(r),phi_est_G(r),phi_est_iS(r)])-0.1 1.5]); % Ensure clear y-axis range
    grid on;
    box off; % Improve aesthetics
    
    % Improve visual clarity with thicker axes
    set(gca, 'LineWidth', 2.5);
    
    % Export settings for publication-quality figures
    figurename=strcat('phi_',sysInfo.name,'N',num2str(sysInfo.N),'d',num2str(sysInfo.d),'M',num2str(obsInfo.M),'L',num2str(length(obsInfo.time_vec)),'sigma',num2str(obsInfo.obs_noise_nsr),'_comparison.png');
    exportgraphics(gcf, figurename, 'BackgroundColor', 'white', 'ContentType', 'vector'); % Save as vector for journal submission

end

%% Visualize \phi (Compatible with SIAM Dynamical Systems Journal & Existing Plot Settings)

% % Define the range for r
% r = linspace(0, rmax, 500); % Smooth interpolation
% 
% % Compute original and estimated phi values
% phi_sys = sysInfo.phi{1}(r);
% phi_L = phi_est_L(r);
% phi_G = phi_est_G(r);
% 
% % % Normalize by Linf norm (maximum absolute value)
% % phi_sys = phi_sys / max(abs(phi_sys));
% % phi_L = phi_L / max(abs(phi_L));
% % phi_B = phi_B / max(abs(phi_B));
% 
% % Normalize by relative l1 norm (mean absolute value)
% % normalize_sys = 1/mean(abs(phi_sys));
% normalize_L = mean(abs(phi_sys(1:200)))/mean(abs(phi_L(1:200)));
% normalize_G = mean(abs(phi_sys(1:200)))/mean(abs(phi_G(1:200)));
% 
% % phi_sys = phi_sys * normalize_sys;
% phi_L = phi_L * normalize_L;
% phi_G = phi_G * normalize_G;
% 
% % Get screen size for figure positioning
% scrsz = [1, 1, 1920, 1080];
% 
% % Create figure with size matching trajectory plots
% phi_fig = figure('Name', 'Comparison of \phi Functions', ...
%                  'NumberTitle', 'off', ...
%                  'Position', [scrsz(3)/8, scrsz(4)/8, scrsz(3)*3/4, scrsz(4)*3/4]);
% 
% % Set background to white
% set(gcf, 'Color', 'w');
% 
% % Define improved color palette
% colors = {[0, 0, 0], [0, 0.447, 0.741], [0.85, 0.325, 0.098]};
% 
% % Plot curves with enhanced styles
% p_true = plot(r, phi_sys, '-', 'LineWidth', 3, 'Color', colors{1}); % True \phi
% hold on;
% p_L = plot(r, phi_L, '--', 'LineWidth', 3, 'Color', colors{2});  % Estimated \phi (L)
% p_G = plot(r, phi_G, '-.', 'LineWidth', 3, 'Color', colors{3});  % Estimated \phi (B)
% % hold off;
% 
% 
% % visualize the results with UQ
% i = L_i;
% ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
% ind = [1:i-1,i+1:nbasis];
% phi_est_LU = @(r) 1/normF(i).*basis_funs{i}(r);
% phi_est_LD = @(r) 1/normF(i).*basis_funs{i}(r);
% for k = 1:nbasis-1
%     if est_L(k) ~= 0
%         phi_est_LU = @(r) phi_est_LU(r) + (est_L(k)+2*errbars_L(k))./normF(ind(k)).*ibasis_funs{k}(r);
%         phi_est_LD = @(r) phi_est_LD(r) + (est_L(k)-2*errbars_L(k))./normF(ind(k)).*ibasis_funs{k}(r);
%     end
% end
% 
% % display the uncertainty region covariance 
% rconf = [0.01:0.01:rmax rmax:-0.01:0.01];
% yconf = [phi_est_LU(0.01:0.01:rmax)*normalize_L phi_est_LD(rmax:-0.01:0.01)*normalize_L];
% 
% % plot(r, phi_est_LU(r)*normalize_L, '-', 'LineWidth', 1, 'Color', colors{2});  % Estimated \phi (L)
% % plot(r, phi_est_LD(r)*normalize_L, '-', 'LineWidth', 1, 'Color', colors{2});  % Estimated \phi (L)
% 
% p = fill(rconf,yconf,colors{2});
% p.EdgeColor = 'none';   
% set(p,'facealpha',.5)
% 
% %
% i = G_i;
% ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
% ind = [1:i-1,i+1:nbasis];
% phi_est_GU = @(r) 1/normF(i).*basis_funs{i}(r);
% phi_est_GD = @(r) 1/normF(i).*basis_funs{i}(r);
% for k = 1:nbasis-1
%     if est_G(k) ~= 0
%         phi_est_GU = @(r) phi_est_GU(r) + (est_G(k)+2*errbars_G(k))./normF(ind(k)).*ibasis_funs{k}(r);
%         phi_est_GD = @(r) phi_est_GD(r) + (est_G(k)-2*errbars_G(k))./normF(ind(k)).*ibasis_funs{k}(r);
%     end
% end
% 
% % display the uncertainty region covariance 
% rconf = [0.01:0.01:rmax rmax:-0.01:0.01];
% yconf = [phi_est_GU(0.01:0.01:rmax).*normalize_G phi_est_GD(rmax:-0.01:0.01).*normalize_G];
% 
% % plot(r, phi_est_BU(r)*normalize_B, '-', 'LineWidth', 1, 'Color', colors{3});  % Estimated \phi (L)
% % plot(r, phi_est_BD(r)*normalize_B, '-', 'LineWidth', 1, 'Color', colors{3});  % Estimated \phi (L)
% 
% p = fill(rconf,yconf,colors{3});
% p.EdgeColor = 'none';   
% set(p,'facealpha',.5)
% 
% % 
% % legend(["True \phi" "Estimated \phi" "Two std region"])
% 
% % plot the density of rho
% yyaxis right                                                                                % display \rho^L_T and its estimator
% axesHandle1         = gca();
% hist1       = histogram(r_rho,100,'EdgeColor',[1 1 1],'Normalization','pdf');
% hist1.FaceAlpha = 0.3;
% 
% hold off;
% 
% % Apply font sizes consistent with your trajectory visualization
% set(gca, 'FontSize', 39, 'FontName', 'Helvetica', 'TickLabelInterpreter', 'latex');
% 
% % Labeling with LaTeX interpreter
% yyaxis left
% xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 50, 'FontName', 'Helvetica');
% ylabel('$\phi(r)$', 'Interpreter', 'latex', 'FontSize', 50, 'FontName', 'Helvetica');
% title(['Estimation of $\phi$ (noise level = ', num2str(noise_to_signal_ratio*100), '\%)'], ...
%     'Interpreter', 'latex', 'FontSize', 50, 'FontName', 'Helvetica');
% 
% % Add a refined legend with SIAM-compatible font size
% legend([p_true, p_L, p_G, hist1],{'True $\phi$', 'Estimated $\hat{\phi}_L$', 'Estimated $\hat{\phi}_G$', '$\rho_r$'}, ...
%     'Interpreter', 'latex', 'FontSize', 33, 'Location', 'northeast');
% 
% % Adjust axes and grid
% xlim([0, rmax]);
% yyaxis left
% ylim([min([phi_sys,phi_est_L(r),phi_est_G(r)])-0.1 max([phi_sys,phi_est_L(r),phi_est_G(r)])+0.1]); % Ensure clear y-axis range
% grid on;
% box off; % Improve aesthetics
% 
% % Improve visual clarity with thicker axes
% set(gca, 'LineWidth', 2.5);
% 
% % Export settings for publication-quality figures
% figurename=strcat('phi_',sysInfo.name,'N',num2str(sysInfo.N),'d',num2str(sysInfo.d),'M',num2str(obsInfo.M),'L',num2str(length(obsInfo.time_vec)),'sigma',num2str(obsInfo.obs_noise_nsr),'_visualization.png');
% exportgraphics(gcf, figurename, 'BackgroundColor', 'white', 'ContentType', 'vector'); % Save as vector for journal submission
% 
% %% visualize coefficients
% switch basis_info.type
%     case 'piecewise'
%         figure
%         for k = 1:nbasis-5
%             line([k,k], [0,coef_L(k)], Color='blue', LineWidth=2)
%             hold on
%             plot(k,coef_L(k),'.',Color = 'blue', MarkerSize=15)
%         end
%         for k = nbasis-4:nbasis
%             line([k-nbasis+20,k-nbasis+20], [0,coef_L(k)], Color='blue')
%             hold on
%             plot(k-nbasis+20,coef_L(k),'.',Color = 'blue',MarkerSize=15)
%         end
% 
%         line([0 20], [0 0], Color='black')
%         hold off
%         title(['FastLaplace, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
%         xticks([1 3 5 7 9 11 15 20])
%         xticklabels({'c_1','c_3','c_5','c_7','c_9','c_{11}','...',strcat('c_{',num2str(nbasis),'}')})
%         ylim([min(coef_L)-0.1 max(coef_L)+0.1])
% %        fontsize(15, "points")
% 
%         figure
%         for k = 1:nbasis-5
%             line([k,k], [0,coef_G(k)], Color='blue', LineWidth=2)
%             hold on
%             plot(k,coef_G(k),'.',Color = 'blue', MarkerSize=15)
%         end
%         for k = nbasis-4:nbasis
%             line([k-nbasis+20,k-nbasis+20], [0,coef_G(k)], Color='blue')
%             hold on
%             plot(k-nbasis+20,coef_G(k),'.',Color = 'blue',MarkerSize=15)
%         end
% 
%         line([0 20], [0 0], Color='black')
%         hold off
%         title(['SparseBayes, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
%         xticks([1 3 5 7 9 11 15 20])
%         xticklabels({'c_1','c_3','c_5','c_7','c_9','c_{11}','...',strcat('c_{',num2str(nbasis),'}')})
%         ylim([min(coef_G)-0.1 max(coef_G)+0.1])
% %        fontsize(15, "points")
% 
% 
%     case {'CSkernel','exp','poly'}
%         figure
%         for k = 1:nbasis
%             line([k,k], [0,coef_L(k)], Color='blue', LineWidth=2)
%             hold on
%             plot(k,coef_L(k),'.',Color = 'blue', MarkerSize=15)
%         end
% 
%         line([0 nbasis+1], [0 0], Color='black')
%         hold off
%         title(['FastLaplace, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
%         xticks([1 3 5 7 9])
%         xticklabels({'c_1','c_3','c_5','c_7','...'})
%         xlim([0 nbasis+1])
%         ylim([min(coef_L)-0.1 max(coef_L)+0.1])
%         %fontsize(15, "points")
% 
%         figure
%         for k = 1:nbasis
%             line([k,k], [0,coef_G(k)], Color='blue', LineWidth=2)
%             hold on
%             plot(k,coef_G(k),'.',Color = 'blue', MarkerSize=15)
%         end
% 
%         line([0 nbasis+1], [0 0], Color='black')
%         hold off
%         title(['SparseBayes, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
%         xticks([1 3 5 7 9])
%         xticklabels({'c_1','c_3','c_5','c_7','...'})
%         xlim([0 nbasis+1])
%         ylim([min(coef_G)-0.1 max(coef_G)+0.1])
% %        fontsize(15, "points")
% 
% end




%% plot components

% plot(0.1:0.01:rmax, 1/normF(i).*basis_funs{i}(0.1:0.01:rmax))
% legend([num2str(i),'-th basis'])
% hold on
% for k = 1:nbasis-1
%     if est(k) ~= 0
%         plot(0.1:0.01:rmax, est(k)./normF(ind(k)).*ibasis_funs{k}(0.1:0.01:rmax),'DisplayName',[num2str(k),'-th basis'])
%     end
% end
% hold off


%% visualize the trajectory

% if sysInfo.d ==1
%     visualize_trajs_1D_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_est);
%     
% elseif sysInfo.d ==2
%         visualize_trajs_2D_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_est);
% 
% end

