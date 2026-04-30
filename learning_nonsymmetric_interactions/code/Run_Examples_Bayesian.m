% Authored by Jinchao Feng and Sui Tang
% RunExamples for the main cases
% clear all;  % commented for batch reproducibility (preserves preset Params)
% close all;

addpaths;

%% Set 2parameters
% if ispc, SAVE_DIR = [getenv('USERPROFILE'), '\DataAnalyses\LearningDynamics']; else, 
%SAVE_DIR = [getenv('HOME'), '/DataAnalyses/LearningDynamics']; 
%end % Please keep this fixed, simply create a symlink ~/DataAnalyses pointing wherever you like                           
VERBOSE                         = 1;                                                                % indicator to print certain output
time_stamp                      = datestr(now, 30);
if ~exist('Params','var'), Params = [];     end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end

%% Load example definitions and let user select one example to run
Examples                        = LoadExampleDefinitions();
ExampleIdx                      = SelectExample(Params, Examples);

%% Get example parameters
Example                         = Examples{ExampleIdx};
sysInfo                        = Example.sysInfo;
solverInfo                     = Example.solverInfo;
obsInfo                        = Example.obsInfo;                                                 % move n to learn_info
obsInfo.VERBOSE                = VERBOSE;
obsInfo.SAVE_DIR               = SAVE_DIR;

%% Learning setup
learnInfo.N = sysInfo.N;
learnInfo.d = sysInfo.d;
learnInfo.order = sysInfo.ode_order;
learnInfo.name = sysInfo.name;
learnInfo.jitter = 1e-6;
learnInfo.dtraj ='false'; % do not compute dtraj during trajectory computation

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

    case 'poly'
        basis_funs = cell(1,nbasis);
        for k = 1:nbasis
            basis_funs{k} = @(x) x.^(k-1);
        end
    
    case 'exp'
        basis_funs = cell(1,nbasis);
        for k = 1:ceil(nbasis)
            basis_funs{k} = @(x) exp(-k.*x);
        end       
end

%%
nT = 1;     %number of trials
errorphis = zeros(2,nT);      %store errors of phis in L-infinity and L2rhoT norms
errortrajs_train = zeros(4,nT);     %store mean and std of trajectory error in training data
errortrajs_test = zeros(4,nT);      %store mean and std of trajectory error in testing data

F_trials = cell(1,nT);

for n = 1:nT
    %% generating data
    [dxpath_test,xpath_test,dxpath_train, xpath_train]=Generate_training_data(sysInfo,obsInfo,solverInfo);
    
    learnInfo.rho_emp = rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);% empirical pairwise distance 
    
    noise_to_signal_ratio = 0;
    obsInfo.obs_noise_nsr = noise_to_signal_ratio; 

    dxpath_true = dxpath_train;

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
    
    %% visualize the training trajectory
    if sysInfo.d ==1
        visualize_training_trajs_1D(sysInfo,obsInfo,solverInfo,learnInfo);

    elseif sysInfo.d ==2
            visualize_training_trajs_2D(sysInfo,obsInfo,solverInfo,learnInfo);

    end


    %% construct F
    F_all = [];
    F_true = [];
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
               elseif strcmp('CSNS',learnInfo.name)
                   for k = 1:nbasis
                       F_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_train(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m)-v_temp,sysInfo.d,[]),2)';
                       F_true_col(k,:) = sum(basis_funs{k}(repmat(r_temp_norm,learnInfo.d,1)).*reshape(dxpath_true(sysInfo.d*sysInfo.N+(i-1)*learnInfo.d+1:sysInfo.d*sysInfo.N+i*learnInfo.d,l,m)-v_temp,sysInfo.d,[]),2)';
                   end
               end
    
               F_all = [F_all F_col];
               F_true = [F_true F_true_col];
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

    %% Sparse Bayesian regression
    % set i-th coefficient to be one
    method1 = 'FastLaplace'; 
    MScriterion1 = ones(1,nbasis)*10^10;

    for i = 1:nbasis
        eta = F(:, i);
        sigma = 0;
        eta = eta + sigma*randn(size(eta,1),size(eta,2));
        rn = randperm(size(F,1));
        Phi_train = -F(rn(1:round(3/4*size(F,1))),[1:i-1,i+1:nbasis]);
        eta_train = eta(rn(1:round(3/4*size(F,1))));
        Phi_test = -F(rn(round(3/4*size(F,1))):end,[1:i-1,i+1:nbasis]);
        eta_test = eta(rn(round(3/4*size(F,1))):end);
        if sum(eta_train~=0) < 0.5*size(eta_train,1)
            continue
        end

        ind = [1:i-1,i+1:nbasis];

        try 
            sigma2 = var(eta)/1e3;
            delta_La = 1e-10;
            lambda_init = [];
            sparsity = [];
            [weights,used_ids,loglik,sigma2,errbars,basis,selected,deleted,alpha,lambda] = FastLaplace(Phi_train,eta_train, sigma2, delta_La,lambda_init,sparsity);


            threshold = [];
            %threshold = 1e-3;
            %threshold = sqrt(sigma2);
            if ~isempty(threshold)
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
    i = L_i;
    eta = F(:, i);
    sigma = 0; % add additional model noise
    eta = eta + sigma*randn(size(eta,1),size(eta,2));
    Phi = -F(:,[1:i-1,i+1:nbasis]);
    ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
    ind = [1:i-1,i+1:nbasis];

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
    phi_est_L = @(r) 1/normF(i).*basis_funs{i}(r);
    for k = 1:nbasis-1
        if est_L(k) ~= 0
            phi_est_L = @(r) phi_est_L(r) + est_L(k)./normF(ind(k)).*ibasis_funs{k}(r);
        end
    end

    coef_L = zeros(nbasis,1);
    coef_L(1:i-1,1) = est_L(1:i-1)./normF(1:i-1)';
    coef_L(i,1) = 1;
    coef_L(i+1:nbasis,1) = est_L(i:nbasis-1)./normF(i+1:nbasis)';

    dr = edges(2) - edges(1);
    phi_fast = @(r) phi_piecewise_const(r, coef_L, rmax);

    result_train = construct_and_compute_traj_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_fast,learnInfo.xpath_train(:,1,:));
    errortrajs_train(:,n) = [result_train.train_traj_error result_train.prediction_traj_error]';
    result_test = construct_and_compute_traj_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_fast,sysInfo.mu0());
    errortrajs_test(:,n) = [result_test.train_traj_error result_test.prediction_traj_error]';

    %% Bayesian Linear Model
    method2 = 'Gaussian'; 
    MScriterion2 = ones(3,nbasis)*10^10;

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

        ind = [1:i-1,i+1:nbasis];

        try 
            [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi_train, eta_train);

            threshold = [];
            if ~isempty(threshold)
                if any(abs(PARAMETER.Value)'.*vecnorm(Phi_train(:,PARAMETER.Relevant))/sqrt(size(Phi,1))<threshold)
                    deleted_ids = (abs(PARAMETER.Value)./norm(PARAMETER.Value)<threshold);
                    PARAMETER.Relevant(deleted_ids) = [];
                    [PARAMETER_t, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi(:,PARAMETER.Relevant), eta);
                    PARAMETER.Relevant = PARAMETER.Relevant(PARAMETER_t.Relevant);
                    PARAMETER.Value = PARAMETER.Value((PARAMETER_t.Relevant));
                end
            end
            
            
            xhat = zeros(length(ind),1);
            xhat(PARAMETER.Relevant) = PARAMETER.Value;
            U = chol(Phi_train(:,PARAMETER.Relevant)'*Phi_train(:,PARAMETER.Relevant)*HYPERPARAMETER.beta + diag(HYPERPARAMETER.Alpha));
            Ui = inv(U);
            SIGMA = Ui * Ui';
            errbars = zeros(length(ind),1);
            errbars(PARAMETER.Relevant) = sqrt(diag(SIGMA));
        catch
            continue
        end

        if sum(Phi_train(:,PARAMETER.Relevant)~=0,"all") < nbasis
            continue
        end

        MScriterion2(1,i) = (sum(errbars(PARAMETER.Relevant).^2)+1./HYPERPARAMETER.beta)./(sum(PARAMETER.Value.^2)+1);
    end
    
    [~,G_i] = min(MScriterion2(1,:)); % find the best 'i-th' basis to set its coefficient to be one
    
    %% BLM with the coefficient of i-th basis set to be one
    i = G_i;
    eta = F(:, i);
    sigma = 0; % add additional model noise
    eta = eta + sigma*randn(size(eta,1),size(eta,2));
    Phi = -F(:,[1:i-1,i+1:nbasis]);
    ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
    ind = [1:i-1,i+1:nbasis];

    [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi, eta);

    threshold = [];
    if ~isempty(threshold)
        if any(abs(PARAMETER.Value)'.*vecnorm(Phi_train(:,PARAMETER.Relevant))/sqrt(size(Phi,1))<threshold)
            deleted_ids = (abs(PARAMETER.Value)./norm(PARAMETER.Value)<threshold);
            PARAMETER.Relevant(deleted_ids) = [];
            [PARAMETER_t, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', Phi(:,PARAMETER.Relevant), eta);
            PARAMETER.Relevant = PARAMETER.Relevant(PARAMETER_t.Relevant);
            PARAMETER.Value = PARAMETER.Value((PARAMETER_t.Relevant));
        end
    end
    
    xhat = zeros(length(ind),1);
    xhat(PARAMETER.Relevant) = PARAMETER.Value;

    U = chol(Phi(:,PARAMETER.Relevant)'*Phi(:,PARAMETER.Relevant)*HYPERPARAMETER.beta + diag(HYPERPARAMETER.Alpha));
    Ui = inv(U);
    SIGMA = Ui * Ui';
    errbars_G = zeros(length(ind),1);
    errbars_G(PARAMETER.Relevant) = sqrt(diag(SIGMA));

    est_G = xhat;
    phi_est_G = @(r) 1/normF(i).*basis_funs{i}(r);
    for k = 1:nbasis-1
        if est_G(k) ~= 0
            phi_est_G = @(r) phi_est_G(r) + est_G(k)./normF(ind(k)).*ibasis_funs{k}(r);
        end
    end

    coef_G = zeros(nbasis,1);
    coef_G(1:i-1,1) = est_G(1:i-1)./normF(1:i-1)';
    coef_G(i,1) = 1;
    coef_G(i+1:nbasis,1) = est_G(i:nbasis-1)./normF(i+1:nbasis)';

    %% ============================================================
    % implicit-SINDy comparison baseline (faithful homogeneous version)
    % Solve F c = 0 by finding the sparsest nonzero vector in an
    % approximate null space using alternating directions method (ADM)
    % ============================================================
    
    run_iSINDy_compare = 1;

    if run_iSINDy_compare
        iS_opts = struct();
        iS_opts.null_dims = 2:nbasis-1;     % candidate null-space dimensions
        iS_opts.lambda_grid = logspace(-4, 0, 20);  % sparsity thresholds for ADM
        iS_opts.maxit = 500;
        iS_opts.tol = 1e-6;
        iS_opts.verbose = false;

        % same train/test split idea as your Bayesian model-selection stage
        rn = randperm(size(F,1));
        train_id = rn(1:round(3/4*size(F,1)));
        test_id  = rn(round(3/4*size(F,1))+1:end);

        F_train = F(train_id, :);
        F_test  = F(test_id,  :);

        [coef_iS, iS_info] = run_implicit_sindy_homogeneous(F_train, F_test, iS_opts);

        % convert coefficient vector into kernel estimator
        phi_est_iS = @(r) 0;
        for k = 1:nbasis
            if abs(coef_iS(k)) > 0
                phi_est_iS = @(r) phi_est_iS(r) + coef_iS(k) ./ normF(k) .* basis_funs{k}(r);
            end
        end

        % fprintf('\nimplicit-SINDy comparison finished.\n');
        % fprintf('  chosen null dimension = %d\n', iS_info.best_null_dim);
        % fprintf('  chosen lambda         = %.4e\n', iS_info.best_lambda);
        % fprintf('  validation residual   = %.4e\n', iS_info.best_residual);
        % fprintf('  nnz(coefficients)     = %d\n', nnz(abs(coef_iS) > 1e-10));
    end
end
%% Visualize \phi
% 
% Define the range for r
r = linspace(0, rmax, 500); % Smooth interpolation

% Compute original and estimated phi values
phi_sys = sysInfo.phi{1}(r);
phi_L = phi_est_L(r);
phi_G = phi_est_G(r);

% Normalization
normalize_L = mean(abs(phi_sys(1:200)))/mean(abs(phi_L(1:200)));
normalize_G = mean(abs(phi_sys(1:200)))/mean(abs(phi_G(1:200)));

% phi_sys = phi_sys * normalize_sys;
phi_L = phi_L * normalize_L;
phi_G = phi_G * normalize_G;
% 

%
errorphis(1,n) = max(abs(phi_sys - phi_L))/max(phi_sys);
errorphis(2,n) = max(abs(phi_sys - phi_G))/max(phi_sys);


% Get screen size for figure positioning
scrsz = [1, 1, 1920, 1080];

% Create figure with size matching trajectory plots
phi_fig = figure('Name', 'Comparison of \phi Functions', ...
                 'NumberTitle', 'off', ...
                 'Position', [scrsz(3)/8, scrsz(4)/8, scrsz(3)*3/4, scrsz(4)*3/4]);

% Set background to white
set(gcf, 'Color', 'w');

% Define improved color palette
colors = {[0, 0, 0], [0, 0.447, 0.741], [0.85, 0.325, 0.098]};

% Plot curves with enhanced styles
phi_true = plot(r, phi_sys, '-', 'LineWidth', 3, 'Color', colors{1}); % True \phi
hold on;
phi_L = plot(r, phi_L, '--', 'LineWidth', 3, 'Color', colors{2});  % Estimated \phi (L)
phi_G = plot(r, phi_G, '-.', 'LineWidth', 3, 'Color', colors{3});  % Estimated \phi (B)
% hold off;


% visualize the results with UQ
i = L_i;
ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
ind = [1:i-1,i+1:nbasis];
phi_est_LU = @(r) 1/normF(i).*basis_funs{i}(r);
phi_est_LD = @(r) 1/normF(i).*basis_funs{i}(r);
for k = 1:nbasis-1
    if est_L(k) ~= 0
        phi_est_LU = @(r) phi_est_LU(r) + (est_L(k)+2*errbars_L(k))./normF(ind(k)).*ibasis_funs{k}(r);
        phi_est_LD = @(r) phi_est_LD(r) + (est_L(k)-2*errbars_L(k))./normF(ind(k)).*ibasis_funs{k}(r);
    end
end

% display the uncertainty region covariance 
rconf = [0.01:0.01:rmax rmax:-0.01:0.01];
yconf = [phi_est_LU(0.01:0.01:rmax)*normalize_L phi_est_LD(rmax:-0.01:0.01)*normalize_L];

p = fill(rconf,yconf,colors{2});
p.EdgeColor = 'none';   
set(p,'facealpha',.5)

%
i = G_i;
ibasis_funs = basis_funs(1,[1:i-1,i+1:nbasis]);
ind = [1:i-1,i+1:nbasis];
phi_est_GU = @(r) 1/normF(i).*basis_funs{i}(r);
phi_est_GD = @(r) 1/normF(i).*basis_funs{i}(r);
for k = 1:nbasis-1
    if est_G(k) ~= 0
        phi_est_GU = @(r) phi_est_GU(r) + (est_G(k)+2*errbars_G(k))./normF(ind(k)).*ibasis_funs{k}(r);
        phi_est_GD = @(r) phi_est_GD(r) + (est_G(k)-2*errbars_G(k))./normF(ind(k)).*ibasis_funs{k}(r);
    end
end

% display the uncertainty region covariance 
rconf = [0.01:0.01:rmax rmax:-0.01:0.01];
yconf = [phi_est_GU(0.01:0.01:rmax).*normalize_G phi_est_GD(rmax:-0.01:0.01).*normalize_G];

p = fill(rconf,yconf,colors{3});
p.EdgeColor = 'none';   
set(p,'facealpha',.5)

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
legend([phi_true, phi_L, phi_G, hist1],{'True $\phi$', 'Estimated $\hat{\phi}_L$', 'Estimated $\hat{\phi}_G$', '$\rho_r$'}, ...
    'Interpreter', 'latex', 'FontSize', 33, 'Location', 'northeast');

% Adjust axes and grid
xlim([0, rmax]);
yyaxis left
ylim([min([phi_sys,phi_est_L(r),phi_est_G(r)])-0.1 max([phi_sys,phi_est_L(r),phi_est_G(r)])+0.1]); % Ensure clear y-axis range
grid on;
box off; % Improve aesthetics

% Improve visual clarity with thicker axes
set(gca, 'LineWidth', 2.5);

% Export settings for publication-quality figures
figurename=strcat('phi_',sysInfo.name,'N',num2str(sysInfo.N),'d',num2str(sysInfo.d),'M',num2str(obsInfo.M),'L',num2str(length(obsInfo.time_vec)),'sigma',num2str(obsInfo.obs_noise_nsr),'_visualization.png');
exportgraphics(gcf, figurename, 'BackgroundColor', 'white', 'ContentType', 'vector'); % Save as vector for journal submission

%% save files
filename=strcat(sysInfo.name,'N',num2str(sysInfo.N),'d',num2str(sysInfo.d),'M',num2str(obsInfo.M),'L',num2str(length(obsInfo.time_vec)),'sigma',num2str(obsInfo.obs_noise_nsr),'.mat');
save(filename, 'sysInfo','obsInfo','learnInfo',"F_trials",'errorphis',"errortrajs_train","errortrajs_test");

%% plot comparison of learned kernel
if run_iSINDy_compare
    % Define the range for r
    r = linspace(0, rmax, 500); % Smooth interpolation

    % Compute original and estimated phi values
    phi_sys = sysInfo.phi{1}(r);
    phi_L = phi_est_L(r);
    phi_G = phi_est_G(r);
    phi_iS = phi_est_iS(r);

    % Normalization
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

%% visualize coefficients
switch basis_info.type
    case 'piecewise'
        figure
        for k = 1:15
            line([k,k], [0,coef_L(k)], Color='blue', LineWidth=2)
            hold on
            plot(k,coef_L(k),'.',Color = 'blue', MarkerSize=15)
        end
        for k = nbasis-4:nbasis
            line([k-nbasis+20,k-nbasis+20], [0,coef_L(k)], Color='blue')
            hold on
            plot(k-nbasis+20,coef_L(k),'.',Color = 'blue',MarkerSize=15)
        end

        line([0 20], [0 0], Color='black')
        hold off
        title(['FastLaplace, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
        xticks([1 3 5 7 9 11 15 100])
        xticklabels({'c_1','c_3','c_5','c_7','c_9','c_{11}','...',strcat('c_{',num2str(nbasis),'}')})
        ylim([min(coef_L)-0.1 max(coef_L)+0.1])
%        fontsize(15, "points")

        figure
        for k = 1:15
            line([k,k], [0,coef_G(k)], Color='blue', LineWidth=2)
            hold on
            plot(k,coef_G(k),'.',Color = 'blue', MarkerSize=15)
        end
        for k = nbasis-4:nbasis
            line([k-nbasis+20,k-nbasis+20], [0,coef_G(k)], Color='blue')
            hold on
            plot(k-nbasis+20,coef_G(k),'.',Color = 'blue',MarkerSize=15)
        end

        line([0 20], [0 0], Color='black')
        hold off
        title(['SparseBayes, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
        xticks([1 3 5 7 9 11 15 20])
        xticklabels({'c_1','c_3','c_5','c_7','c_9','c_{11}','...',strcat('c_{',num2str(nbasis),'}')})
        ylim([min(coef_G)-0.1 max(coef_G)+0.1])
%        fontsize(15, "points")


    case {'exp','poly'}
        figure
        for k = 1:nbasis
            line([k,k], [0,coef_L(k)], Color='blue', LineWidth=2)
            hold on
            plot(k,coef_L(k),'.',Color = 'blue', MarkerSize=15)
        end

        line([0 nbasis+1], [0 0], Color='black')
        hold off
        title(['FastLaplace, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
        xticks([1 3 5 7 9])
        xticklabels({'c_1','c_3','c_5','c_7','...'})
        xlim([0 nbasis+1])
        ylim([min(coef_L)-0.1 max(coef_L)+0.1])
        %fontsize(15, "points")

        figure
        for k = 1:nbasis
            line([k,k], [0,coef_G(k)], Color='blue', LineWidth=2)
            hold on
            plot(k,coef_G(k),'.',Color = 'blue', MarkerSize=15)
        end

        line([0 nbasis+1], [0 0], Color='black')
        hold off
        title(['SparseBayes, \nsr = ',num2str(noise_to_signal_ratio)],'FontSize',12)
        xticks([1 3 5 7 9])
        xticklabels({'c_1','c_3','c_5','c_7','...'})
        xlim([0 nbasis+1])
        ylim([min(coef_G)-0.1 max(coef_G)+0.1])
end

%% visualize the trajectory

if sysInfo.d ==1
    visualize_trajs_1D_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_fast);

elseif sysInfo.d ==2
        visualize_trajs_2D_NS(sysInfo,obsInfo,solverInfo,learnInfo,phi_fast);

end




