clear all;
addpaths;
load("FM5L3M3");
%load("FM10L3M3");
%load("FM12L3M3");
%load("FM14L3M3");
%load("FM16L3M3");
%load("FM18L3M3");
%load("FM20L3M3");
%load("FM22L3M3");
%load("AD10L3M4.mat");
%load("FM10nu52");


learnInfo.v = 5/2;

learnInfo.hyp = [1, 1, 1, 1, 0.5, 1, 1];

M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;

%Adjustable parameters.
CGErrorTol = 10^(-5);
CG_ITER_LIMIT = 100;
mFORGLIK = 100;
lFORGLIK = 100;
rangeOfI = 4;
jitter = 10^(-6);
HVAL = 10^(-5);
GlikRuns = 30;


%Decomp for K_E.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train(D*n+1:2*D*n,:,:);
%dataA = learnInfo.dxpath_train(1:D*n,:,:);
% dataAAA = learnInfo.dxpath_train(D*n+1:2*D*n,:,:);




nT = 3;     %number of trials
errorphis = zeros(4,nT);      %store errors of phis in L-infinity and L2rhoT norms
errortrajs_train = zeros(4,nT);     %store mean and std of trajectory error in training data
errortrajs_test = zeros(4,nT);      %store mean and std of trajectory error in testing data
hypparameters = zeros(length(learnInfo.hyp0),nT);   %store estimated hyperparameters
errorhyp = zeros(3,nT);      %store hyperparameter errors

runtimes = zeros(4,nT);





for k = 1 : nT
    %Get constants.
    
    deltaE = exp(learnInfo.hyp(1));
    omegaE = exp(learnInfo.hyp(2));
    deltaA = exp(learnInfo.hyp(3));
    omegaA = exp(learnInfo.hyp(4));
    sigma = exp(learnInfo.hyp(5))^2;
    
    %If sigma is NaN, there is no noise. Use jitter factor.
    if isnan(sigma)
        sigma = 10^(-6);
    end
    
    M = obsInfo.M;
    LForDecomp = obsInfo.L;
    n = learnInfo.N;
    D = learnInfo.d;
    
    
    
    
    tic
    learnInfo.option = 'alldata';
    [TRUEfval, TRUEdfval,~] = Glik(learnInfo,learnInfo.hyp)
    toc
    
    
    
    tic
    [fval2, dfval2,~] = GlikSTEPreConSTE52(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, @IdentityPrecon)
    toc
    
    
    
    
    %Full Glik
    
    
    tic
    learnInfo.hyp = [1, 1, 1, 1, 0.5, 1, 1];
    learnInfo.option = 'alldata';
    [TRUEfval, TRUEdfval,~] = Glik(learnInfo,learnInfo.hyp)
    Glik_hyp = @(hyp)Glik(learnInfo,hyp);
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -GlikRuns);
    
    toc
    
    tic
    learnInfo.option = 'alldata';
    [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);
    learnInfo.invK = pinv(learnInfo.K);
    learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;
    
    
    hypparameters(:,k) = exp(learnInfo.hyp)
    if strcmp('ODS',learnInfo.name) 
        hypparameters(5:end,k) = learnInfo.hyp(5:end);
    end
    
    learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'E');
    learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'A');
    
    
    toc
    
    
    
    tic
    range = [0, learnInfo.rhoLT.edges(max(find(learnInfo.rhoLT.rdens~=0)))];
        
    [learnInfo, errorphis(1,k),errorphis(2,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'E')
    [learnInfo, errorphis(3,k),errorphis(4,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'A')
    result_train = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, learnInfo.xpath_train(:,1,:));
    errortrajs_train(:,k) = [result_train.train_traj_error result_train.prediction_traj_error]'
    result_test = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo,sysInfo.mu0());
    errortrajs_test(:,k) = [result_test.train_traj_error result_test.prediction_traj_error]'
    toc
    
    
    
    
    
    %START APPROX
    
    
    tic
    
    learnInfo.hyp = [1, 1, 1, 1, 0.5, 1, 1];
    [fval, dfval,~] = GlikSTEPreConSTE52(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, @IdentityPrecon);
    
    Glik_hyp = @(hyp)GlikSTEPreConSTE52(learnInfo, hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, @IdentityPrecon);
    
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -GlikRuns);
    
    
    %If accelerated, then we only need Ym.
    learnInfo.option = 'alldata';
    learnInfo = GetYm(learnInfo,learnInfo.hyp);
    
        
    hypparameters(:,k) = exp(learnInfo.hyp)
    if strcmp('ODS',learnInfo.name) 
        hypparameters(5:end,k) = learnInfo.hyp(5:end);
    end
    
    deltaE = exp(learnInfo.hyp(1));
    omegaE = exp(learnInfo.hyp(2));
    deltaA = exp(learnInfo.hyp(3));
    omegaA = exp(learnInfo.hyp(4));
    sigma = exp(learnInfo.hyp(5))^2;
    
    %If sigma is NaN, there is no noise. Use jitter factor.
    if isnan(sigma)
        sigma = 10^(-6);
    end
    
    M = obsInfo.M;
    LForDecomp = obsInfo.L;
    n = learnInfo.N;
    D = learnInfo.d;
    
    
    %TODO: Add Preconditioner.
    PreConInv = @(x) x;
    % [~, PreConInvMat, ~] = ConstructNystPreconNoK(learnInfo, LForDecomp, M, rangeOfI, jitter);
    % PreConInv = @(x) PreConInvMat * x;
    
    
    
    
    [U_sE, R_sE, ~, ~] = TotalDecompForDebug52(data, data, omegaE, deltaE, n, D, M, LForDecomp);
    
    KE = U_sE * R_sE * U_sE';
    
    MultByKE = @(x) KE * x;
    
    [U_sA, R_sA, ~, ~] = TotalDecompForDebug52(data, dataA, omegaA, deltaA, n, D, M, LForDecomp);
    
    KA = U_sA * R_sA * U_sA';
    
    MultByKA = @(x) KA * x;
    
    
    
    %Multiply by K + sigmaI.
    MultByWholeK = @(x) MultByKE(x) + MultByKA(x) + sigma*x;
    
    %Multiply by entire (K+sigmaI)^-1.
    multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, lFORGLIK);
    
    
    %Now let's test matrix multiplication.
    multMatByKInv = @(X) RunCGOnMatrixInitGuesser(MultByWholeK, X, PreConInv, CGErrorTol, CG_ITER_LIMIT);
    
    learnInfo.invKTimesYm = multByKInv(learnInfo.Ym);
    toc
    
    
    tic
    visualize_phis_CG(sysInfo,obsInfo,learnInfo,'E', multMatByKInv);
    visualize_phis_CG(sysInfo,obsInfo,learnInfo,'A', multMatByKInv);
    toc
    
    
    tic
    range = [0, learnInfo.rhoLT.edges(max(find(learnInfo.rhoLT.rdens~=0)))];
        
    [learnInfo, errorphis(1,k),errorphis(2,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'E')
    [learnInfo, errorphis(3,k),errorphis(4,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'A')
    result_train = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, learnInfo.xpath_train(:,1,:));
    errortrajs_train(:,k) = [result_train.train_traj_error result_train.prediction_traj_error]'
    result_test = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo,sysInfo.mu0());
    errortrajs_test(:,k) = [result_test.train_traj_error result_test.prediction_traj_error]'
    toc



    


