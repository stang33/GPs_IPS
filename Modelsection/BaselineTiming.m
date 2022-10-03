clear all;
addpaths;
%load("FM5L3M3");
%load("FM10L3M3");
%load("FM12L3M3");
%load("FM14L3M3");
%load("FM16L3M3");
%load("FM18L3M3");
%load("FM20L3M3");
%load("FM22L3M3");
%load("AD10L3M4.mat");
%load("FM10nu52");
%load("FM40L3M2");

load("FM40v2.mat");


learnInfo.v = 5/2;

learnInfo.hyp = [1, 1, 1, 1, 0.5, 1, 1];

M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;

%Adjustable parameters.
CGErrorTol = 10^(-6);
CG_ITER_LIMIT = 150;
mFORGLIK = 150;
lFORGLIK = 150;
rangeOfI = 4;
jitter = 10^(-6);
HVAL = 10^(-5);
GlikRuns = 100;


%Decomp for K_E.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train(D*n+1:2*D*n,:,:);
%dataA = learnInfo.dxpath_train(1:D*n,:,:);
% dataAAA = learnInfo.dxpath_train(D*n+1:2*D*n,:,:);




nT = 2;     %number of trials
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
    
    
    
    
    tic;
    learnInfo.option = 'alldata';
    [TRUEfval, TRUEdfval,~] = Glik(learnInfo,learnInfo.hyp);
    runtimes(1,k) = toc;
    
    
    
    %Full Glik
    
    learnInfo.hyp = [1, 1, 1, 1, 0.5, 1, 1];
    learnInfo.option = 'alldata';
    [TRUEfval, TRUEdfval,~] = Glik(learnInfo,learnInfo.hyp)
    Glik_hyp = @(hyp)Glik(learnInfo,hyp);
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -GlikRuns);
    
    runtimes(2,k) = toc;
    


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
    runtimes(3,k) = toc;
    
    
    

    range = [0, learnInfo.rhoLT.edges(max(find(learnInfo.rhoLT.rdens~=0)))];
        
    [learnInfo, errorphis(1,k),errorphis(2,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'E');
    [learnInfo, errorphis(3,k),errorphis(4,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'A');
    result_train = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, learnInfo.xpath_train(:,1,:));
    errortrajs_train(:,k) = [result_train.train_traj_error result_train.prediction_traj_error]';
    result_test = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo,sysInfo.mu0());
    errortrajs_test(:,k) = [result_test.train_traj_error result_test.prediction_traj_error]';
    runtimes(4,k) = toc;

    save("postFullIter");

end




avgtime = zeros(4,1);
stdtime = zeros(4,1);

for i = 1 : 4
    avgtime(i,1) = mean(runtimes(i,:));
    stdtime(i,1) = std(runtimes(i,:));
end

avgtime
stdtime


avgkerror = zeros(4,1);
stdkerror = zeros(4,1);
avgetrain = zeros(4,1);
stdetrain = zeros(4,1);
avgetest = zeros(4,1);
stdetest = zeros(4,1);

for i = 1 : 4
    avgkerror(i,1) = mean(errorphis(i,:));
    stdkerror(i,1) = std(errorphis(i,:));
    avgetrain(i,1) = mean(errortrajs_train(i,:));
    stdetrain(i,1) = std(errortrajs_train(i,:));
    avgetest(i,1) = mean(errortrajs_test(i,:));
    stdetest(i,1) = std(errortrajs_test(i,:));
end

avgkerror
stdkerror
avgetrain
stdetrain
avgetest
stdetest


save("postFull");
