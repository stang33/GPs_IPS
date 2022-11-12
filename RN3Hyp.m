function RN3Hyp(MVAL)
addpaths;

%Datagen stuff.

fullname = strcat("FM_RNG_M",num2str(MVAL));

if isfile(strcat(fullname,".mat"))
     % File exists.
     load(fullname);
     fprintf('\n Using previously generated data ......');
else
     % File does not exist.

     fprintf('\n No such file. Generating new data ......');

sysInfo                        = FM_def(MVAL).sysInfo;
solverInfo                     = FM_def(MVAL).solverInfo;
obsInfo                        = FM_def(MVAL).obsInfo;                                                 % move n to learn_info
obsInfo.MrhoT = 2000;
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

learnInfo.v = 5/2; % 1/2, 3/2, 5/2, 7/2

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


%Now, data has been loaded. Run the standard RN stuff.


%Enter method to use for greatest likelihood.
%@GlikSTEPreConSTE52 - the full accelerated method.
%@GlikSTEPreConSTE52EK - Use exact derivative calculations.
%@GlikSTEPreConSTE52EFVAL - Use exact fval calculation.
GlikMethod = @GlikNo1to4;

%Preconditioner method.
PreconMethod = @RandomNystbackup2;

%Constants.
learnInfo.v = 5/2;
M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;

%Adjustable parameters.
CGErrorTol = 10^(-6);   %Error tolerance in the CG algorithm.
CG_ITER_LIMIT = 150;   %Max number of CG iterations.
mFORGLIK = CG_ITER_LIMIT;   %Number of test vectors for Lanczos.
lFORGLIK = CG_ITER_LIMIT;   %Number of CG iters from Lanczos.
rangeOfI = floor(M*LForDecomp*n*D / 30);   %Size of preconditioner.
jitter = 10^(-5);   %Value to use for sigma if sigma is too small.
HVAL = 10^(5);   %NOT USED.
GlikRuns = 100;   %Runs of minimizing.
alpha = 0.4;  %Range of randomness in intialization of hyps.
nT = 10;     %Number of trials.

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
hyp_errors_rel = zeros(3,nT);      %store hyperparameter errors
hyp_errors_abs = zeros(3,nT);      %store hyperparameter errors


for k = 1 : nT

    %These are guesses for hyperparameters so we can center
    %intitialization values for the learned hyperparameters.
    originalHyps = [1, 1, 1, 1, .01, 1, 1];

    %Add randomness.
    for i = 1 : 7
        originalHyps(i) = originalHyps(i) + alpha * (2 * rand(1,1) - 1);
    end
    learnInfo.hyp = originalHyps;

    %Get constants from hyperparameters.
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
    [fval2, dfval2,~] = GlikMethod(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, PreconMethod)
    runtimes(1,k) = toc;

    
    %Reset and run the actual optimization.
    learnInfo.hyp = originalHyps;
    [fval, dfval,~] = GlikMethod(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, PreconMethod);
    Glik_hyp = @(hyp)GlikMethod(learnInfo, hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI, PreconMethod);
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -GlikRuns);
    runtimes(2,k) = toc;

    learnInfo.hyp(1:4) = [1,1,1,1];

    hyp_errors_rel(1,k) = abs(exp(learnInfo.hyp(5))^2 - .01) / abs(.01);
    hyp_errors_rel(2,k) = abs(exp(learnInfo.hyp(6)) - 1.5) / abs(1.5);
    hyp_errors_rel(3,k) = abs(exp(learnInfo.hyp(7)) - 0.5) / abs(0.5);

    hyp_errors_abs(1,k) = abs(exp(learnInfo.hyp(5))^2 - .01);
    hyp_errors_abs(2,k) = abs(exp(learnInfo.hyp(6)) - 1.5);
    hyp_errors_abs(3,k) = abs(exp(learnInfo.hyp(7)) - 0.5);


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
    


    learnInfo.option = 'alldata';
    [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);
    learnInfo.invK = pinv(learnInfo.K);
    learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;
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

    %Save file.
    filename = strcat("RN3HypITER",num2str(MVAL));
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

hyp_errors_rel
hyp_errors_abs


%Save final report.
filename = strcat("RN3HypRESULT",num2str(MVAL));
save(filename);





end


