function Run_prediction_Sui(MVAL)
%Usage: Call RNystRun() and enter the desired M-value as a parameter.
%The data included will be automatically loaded to specification.
%% this is for exact hyperameter and only M varing 

addpaths;


fullname = strcat("FMN10_RNG_M",num2str(MVAL));

if isfile(strcat(fullname,".mat"))
     % File exists.
     load(fullname);
     fprintf('\n Using previously generated data ......');
else
     % File does not exist.

     fprintf('\n No such file. Generating new data ......');

sysInfo                        = FM_def().sysInfo;
solverInfo                     = FM_def().solverInfo;
obsInfo                        = FM_def().obsInfo;                                                 % move n to learn_info
obsInfo.MrhoT = 1;
saveON= 0;
plotON = 0;
learnInfo.rhoLT = Generate_rhoT(sysInfo,obsInfo,solverInfo,saveON,plotON);% empirical pairwise distance

if obsInfo.obs_noise>0
  obsInfo.use_derivative     = true;
end

learnInfo.N = sysInfo.N;
learnInfo.d = sysInfo.d;
learnInfo.order = sysInfo.ode_order;
learnInfo.name = sysInfo.name;
learnInfo.jitter = 1e-6;
learnInfo.Cov = 'Matern';
learnInfo.dtraj ='false'; % do not compute dtraj during trajectory computation

learnInfo.v = 3/2; % 1/2, 3/2, 5/2, 7/2

learnInfo.hyp0 = log([rand(1,7)]);

[dxpath_test,xpath_test,dxpath_train, xpath_train]=Generate_training_data(sysInfo,obsInfo,solverInfo);
learnInfo.rho_emp = rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);% empirical pairwise distance 
dxpath_train= trajUnifNoiseAdditive(dxpath_train, obsInfo.obs_noise);



X_train=[];
Y_train=[];


for i = 1:obsInfo.M
    for j = 1:size(xpath_train,2)
        X_train = [X_train;xpath_train(:,j,i)];
        if sysInfo.ode_order ==2
            Y_train = [Y_train;dxpath_train(sysInfo.d*sysInfo.N+1:end,j,i)];
        else
            Y_train = [Y_train;dxpath_train(:,j,i)];
        end
    end
end

learnInfo.L= size(xpath_train,2)*size(xpath_train,3);
learnInfo.X = X_train;
learnInfo.Y = Y_train;

learnInfo.xpath_train = xpath_train;
learnInfo.dxpath_train = dxpath_train;

learnInfo.hyp = learnInfo.hyp0;   

save(fullname);

fprintf('\n Done generating! Begin to learn ......');

end

%Enter method to use for greatest likelihood.
%@GlikSTEPreConSTE52 - the full accelerated method.
%@GlikSTEPreConSTE52EK - Use exact derivative calculations.
%@GlikSTEPreConSTE52EFVAL - Use exact fval calculation.
GlikMethod = @GlikNo1to4;

%Preconditioner method.
PreconMethod = @RandomNystbackup2_Sui;

%Load data. Change this as needed to call your exact dataset.
%firstpart = strcat("FM1011","5");
%fullname = strcat(firstpart, "FixL.mat");
% fullname = "FM600";
% fullname='FM10113FixL.mat';
% load(fullname);

%Constants.
learnInfo.v = 3/2;
M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;

%Adjustable parameters.
CGErrorTol = 10^(-12);   %Error tolerance in the CG algorithm.
CG_ITER_LIMIT = 100;   %Max number of CG iterations.
mFORGLIK = CG_ITER_LIMIT;   %Number of test vectors for Lanczos.
lFORGLIK = CG_ITER_LIMIT;   %Number of CG iters from Lanczos.
rangeOfI = 40;   %Size of preconditioner.
jitter = 10^(-4);   %Value to use for sigma if sigma is too small.
HVAL = 10^(-5);   %NOT USED.
nT = 1;     %Number of trials.

%Decomp for decompositon.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train(D*n+1:2*D*n,:,:);

%Result storage containers.
errorphis = zeros(4,nT);      %store errors of phis in L-infinity and L2rhoT norms
errortrajs_train = zeros(4,nT);     %store mean and std of trajectory error in training data
errortrajs_test = zeros(4,nT);      %store mean and std of trajectory error in testing data
hypparameters = zeros(length(learnInfo.hyp0),nT);   %store estimated hyperparameters
errorhyp = zeros(3,nT);      %store hyperparameter errors
runtimes = zeros(4,nT);      %store runtimes of each section


for k = 1 : nT

    %These are guesses for hyperparameters so we can center
    %intitialization values for the learned hyperparameters.
    originalHyps = [25./log(obsInfo.M), 1.6, 23./log(obsInfo.M), 1.2, .1, 1.5, 0.5];% noise 0.1^2=0.01

    %Add randomness.
    
        originalHyps(5) = originalHyps(5) + 0.0024 * (2 * rand(1,1) - 1);
        originalHyps(6) = originalHyps(6) + 0.05 * (2 * rand(1,1) - 1);
        originalHyps(7) = originalHyps(7) + 0.001 * (2 * rand(1,1) - 1);


    learnInfo.hyp = log(originalHyps);



    %Time a single run of Greatest Likelihood.
%     tic;
%     [fval2, dfval2,~] = GlikMethod(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, PreconMethod)
%     runtimes(1,k) = toc;
% 
%     
%     %Reset and run the actual optimization.
%     learnInfo.hyp = originalHyps;
%     [fval, dfval,~] = GlikMethod(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, PreconMethod);
%     Glik_hyp = @(hyp)GlikMethod(learnInfo, hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, PreconMethod);
%     [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -GlikRuns);
%     runtimes(2,k) = toc;

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
    
    %Now we set up for prediction.
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
    
    %One more preconditioner call.
    [~, PreConInvRaw] = PreconMethod(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma);
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
%     result_train = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, learnInfo.xpath_train(:,1,:));
%     errortrajs_train(:,k) = [result_train.train_traj_error result_train.prediction_traj_error]';
%     result_test = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo,sysInfo.mu0());
%     errortrajs_test(:,k) = [result_test.train_traj_error result_test.prediction_traj_error]';
    runtimes(4,k) = toc;

    %Save file.
    filename = strcat("RNystITERE",num2str(MVAL));
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


%Save final report.
filename = strcat("RNystResultE",num2str(MVAL));
save(filename);


end
