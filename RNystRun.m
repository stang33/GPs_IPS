function RNystRun(MVAL)
%Usage: Call RNystRun() and enter the desired M-value as a parameter.
%The data included will be automatically loaded to specification.

addpaths;

%Load data.
firstpart = strcat("FM1011",num2str(MVAL));
fullname = strcat(firstpart, "FixL.mat");
load(fullname);

%Constants.
learnInfo.v = 5/2;
M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;

%Adjustable parameters.
CGErrorTol = 10^(-6);
CG_ITER_LIMIT = 150;
mFORGLIK = CG_ITER_LIMIT;
lFORGLIK = CG_ITER_LIMIT;
rangeOfI = floor(M*LForDecomp*n*D / 10);
jitter = 10^(-5);
HVAL = 10^(-5);
GlikRuns = 100;
alpha = 1;
nT = 4;     %number of trials



%Decomp for K_E.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train(D*n+1:2*D*n,:,:);

errorphis = zeros(4,nT);      %store errors of phis in L-infinity and L2rhoT norms
errortrajs_train = zeros(4,nT);     %store mean and std of trajectory error in training data
errortrajs_test = zeros(4,nT);      %store mean and std of trajectory error in testing data
hypparameters = zeros(length(learnInfo.hyp0),nT);   %store estimated hyperparameters
errorhyp = zeros(3,nT);      %store hyperparameter errors
runtimes = zeros(4,nT);


for k = 1 : nT


%     originalHyps = [1, 1, 1, 1, 0.5, 1, 1];
%     originalHyps = 2.5 + alpha * (2 * rand(1,7) - 1);
%     originalHyps(3) = 1.6 + alpha * (2 * rand(1,1) - 1);
%     originalHyps(4) = 2.3 + alpha * (2 * rand(1,1) - 1);
%     originalHyps(5) = .01 + alpha * (2 * rand(1,1) - 1);
%     originalHyps(6) = 1.5 + alpha * (2 * rand(1,1) - 1);
%     originalHyps(7) = 0.5 + alpha * (2 * rand(1,1) - 1);

    %We use the theoretically harder starting parameters, as they seem
    % to have a lower failure rate.
    originalHyps = 1 + alpha * (2 * rand(1,7) - 1);    
    learnInfo.hyp = originalHyps;




    %Get constants.
    deltaE = exp(learnInfo.hyp(1));
    omegaE = exp(learnInfo.hyp(2));
    deltaA = exp(learnInfo.hyp(3));
    omegaA = exp(learnInfo.hyp(4));
    sigma = exp(learnInfo.hyp(5))^2;
    
    %If sigma is NaN, there is no noise. Use jitter factor.
    if isnan(sigma)
        sigma = jitter;
    end

    %Time a single run of Greatest Likelihood.
    tic;
    [fval2, dfval2,~] = GlikSTEPreConSTE52(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, @RandomNyst)
    runtimes(1,k) = toc;

    
    %Reset and run the actual optimization.
    learnInfo.hyp = originalHyps;
    [fval, dfval,~] = GlikSTEPreConSTE52(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, @RandomNyst);
    Glik_hyp = @(hyp)GlikSTEPreConSTE52(learnInfo, hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, @RandomNyst);
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -GlikRuns);
    runtimes(2,k) = toc;

    
    %If accelerated, then we only need Ym to predict.
    learnInfo.option = 'alldata';
    learnInfo = GetYm(learnInfo,learnInfo.hyp);
    
    %Set final hyperparameters.   
    hypparameters(:,k) = exp(learnInfo.hyp);
    deltaE = exp(learnInfo.hyp(1));
    omegaE = exp(learnInfo.hyp(2));
    deltaA = exp(learnInfo.hyp(3));
    omegaA = exp(learnInfo.hyp(4));
    sigma = exp(learnInfo.hyp(5))^2;
    
    %If sigma is NaN, there is no noise. Still use jitter factor.
    if isnan(sigma)
        sigma = jitter;
    end
    
    
    X = learnInfo.X;
    dN = learnInfo.d*learnInfo.N*learnInfo.order;
    L = length(X)/dN;
    LForDecomp = L / M;
    
    %Decompose for the final kernel methods for prediction.
    [U_sE, R_sE, ~, ~] = TotalDecompForDebug52(data, data, omegaE, deltaE, n, D, M, LForDecomp); 
    KE = U_sE * R_sE * U_sE';
    MultByKE = @(x) KE * x;
    
    [U_sA, R_sA, ~, ~] = TotalDecompForDebug52(data, dataA, omegaA, deltaA, n, D, M, LForDecomp);
    KA = U_sA * R_sA * U_sA';
    MultByKA = @(x) KA * x;
    
    [~, PreConInvRaw] = RandomNyst(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma);
    PreConInv = @(x) PreConInvRaw * x;

    %Multiply by K + sigmaI.
    MultByWholeK = @(x) MultByKE(x) + MultByKA(x) + sigma*x;
    
    %Multiply by entire (K+sigmaI)^-1.
    multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, lFORGLIK);
    
    %Matrix multiplication.
    multMatByKInv = @(X) RunCGOnMatrixInitGuesser(MultByWholeK, X, PreConInv, CGErrorTol, CG_ITER_LIMIT);
    
    learnInfo.invKTimesYm = multByKInv(learnInfo.Ym);

    %Visualize kernel.
    visualize_phis_CG(sysInfo,obsInfo,learnInfo,'E', multMatByKInv);
    visualize_phis_CG(sysInfo,obsInfo,learnInfo,'A', multMatByKInv);
    runtimes(3,k) = toc;
    
    
    %Calculate errors.
    range = [0, learnInfo.rhoLT.edges(max(find(learnInfo.rhoLT.rdens~=0)))];
    [learnInfo, errorphis(1,k),errorphis(2,k)] = errornorms_phis_CG(sysInfo,obsInfo,learnInfo,range,'E', multMatByKInv);
    [learnInfo, errorphis(3,k),errorphis(4,k)] = errornorms_phis_CG(sysInfo,obsInfo,learnInfo,range,'A', multMatByKInv);
    result_train = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, learnInfo.xpath_train(:,1,:));
    errortrajs_train(:,k) = [result_train.train_traj_error result_train.prediction_traj_error]';
    result_test = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo,sysInfo.mu0());
    errortrajs_test(:,k) = [result_test.train_traj_error result_test.prediction_traj_error]';
    runtimes(4,k) = toc;

    filename = strcat("RNystITER",num2str(MVAL));
    save(filename);

end

%Print results for ease of use.
runtimes


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

hypparameters



filename = strcat("RNystResult",num2str(MVAL));

save(filename);


end
