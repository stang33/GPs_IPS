clear all;
addpaths;
load("FM10L3M3");
%load("FM20L3M3");
%load("AD10L3M4.mat");
%load("FM10nu52");

M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;

%Adjustable parameters.
CGErrorTol = 10^(-5);
CG_ITER_LIMIT = 100;
mFORGLIK = 50;
lFORGLIK = 50;
rangeOfI = 4;
jitter = 10^(-6);
HVAL = 10^(-5);


%Decomp for K_E.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train(D*n+1:2*D*n,:,:);
%dataA = learnInfo.dxpath_train(1:D*n,:,:);
% dataAAA = learnInfo.dxpath_train(D*n+1:2*D*n,:,:);

USE_ACCELERATION = true;




if ~USE_ACCELERATION

learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'E');
learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'A');

else

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










% trueK = learnInfo.K;
% 
% [U_s, R_s, ~, ~] = TotalDecompForDebug(data, omegaE, deltaE, n, D, M, LForDecomp, learnInfo);
% [U_sA, R_sA, ~, ~] = TotalDecompForDebugA(data, dataA, omegaA, deltaA, n, D, M, LForDecomp);
% 
% out = U_s*R_s*U_s' - trueK;
% %out(1:20,1:20)
% 
% out2 = U_sA*R_sA*U_sA';
% %out2(1:20,1:20)
% 
% fullOut = U_s*R_s*U_s' + U_sA*R_sA*U_sA' + sigma*eye(n*D*M*LForDecomp);
% 
% 
% 
% %This norm is close.
% norm(fullOut - trueK)
% 
% 



%TODO: Add Preconditioner.
%PreConInv = @(x) x;
% [logDetPreCon, PreConInvMat, PreCon] = ConstructNystPreconNoK(learnInfo, LForDecomp, M, rangeOfI, jitter);
% PreConInv = @(x) PreConInvMat * x;


% 
% %TODO: For now, ID Precon.
PreConInv = @(x) x;
logDetPreCon = 0;



[U_re, P_r, P_c, rhoVect, ~, newRawIndices, ~, ms, ls] = SparseDecomp(data, data, learnInfo, M, LForDecomp, omegaE, deltaE);




%Multiply by K_E.
MultByKNoNoiseTerm = @(x)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_re, P_r, P_c, rhoVect, newRawIndices, deltaE, ls, ms);

%Decomp for K_A.
[U_re, P_r, P_c, rhoVect, ~, newRawIndices, ds, ms, ls] = SparseDecomp(data, dataA, learnInfo, M, LForDecomp, omegaA, deltaA);

%Multiply by K_A.
MultByKNoNoiseTermA = @(x)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_re, P_r, P_c, rhoVect, newRawIndices, deltaA, ls, ms);

%Multiply by K + sigmaI.
MultByWholeK = @(x) MultByKNoNoiseTerm(x) + MultByKNoNoiseTermA(x) + sigma*x;

%Multiply by entire (K+sigmaI)^-1.
multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, CG_ITER_LIMIT);

%Now let's test matrix multiplication.
multMatByKInv = @(X) RunCGOnMatrixInitGuesser(MultByWholeK, X, PreConInv, CGErrorTol, CG_ITER_LIMIT);



% 
% %Compute the Ym product once.
% learnInfo.invKTimesYm = multByKInv(learnInfo.Ym);











% trueK = learnInfo.K;
% origCon = cond(trueK)
% newCon = cond(PreConInv(trueK))
% norm(trueK - PreCon)
% norm(pinv(trueK) - PreConInvMat)
% norm(pinv(trueK) - pinv(PreCon))
% norm(PreConInvMat - pinv(PreCon))




tic
learnInfo.option = 'alldata';
[TRUEfval, TRUEdfval,~] = Glik(learnInfo,learnInfo.hyp)
    Glik_hyp = @(hyp)Glik(learnInfo,hyp);
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -40);
    

    %If not accelerated, we must run additonal Glik to make full K and Ym.
    if ~USE_ACCELERATION
        learnInfo.option = 'alldata';
        [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);
        learnInfo.invK = pinv(learnInfo.K);
        learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;
    else
        %If accelerated, then we only need Ym.
        learnInfo.option = 'alldata';
        learnInfo = GetYm(learnInfo,learnInfo.hyp);
    end

    
    hypparameters(:,k) = exp(learnInfo.hyp)
    if strcmp('ODS',learnInfo.name) 
        hypparameters(5:end,k) = learnInfo.hyp(5:end);
    end

    learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'E');
learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'A');


toc

% tic
% learnInfo.option = 'subset';  % doesn't work for ODS
% learnInfo.Nsub = 8;
% learnInfo.sub = randsample(1:learnInfo.N,learnInfo.Nsub);
% [TRUEfval, TRUEdfval,~] = Glik(learnInfo,learnInfo.hyp)
% toc
% 
% 
% tic
% [fval, dfval,~] = GlikSTEPreConSTE(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI)
% toc

tic
[fval, dfval,~] = GlikNOPRETEST(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI)

 Glik_hyp = @(hyp)GlikNOPRETEST(learnInfo, hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI);

    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -40);
    


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

    PreConInv = @(x) x;
logDetPreCon = 0;



[U_re, P_r, P_c, rhoVect, ~, newRawIndices, ~, ms, ls] = SparseDecomp(data, data, learnInfo, M, LForDecomp, omegaE, deltaE);




%Multiply by K_E.
MultByKNoNoiseTerm = @(x)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_re, P_r, P_c, rhoVect, newRawIndices, deltaE, ls, ms);

%Decomp for K_A.
[U_re, P_r, P_c, rhoVect, ~, newRawIndices, ds, ms, ls] = SparseDecomp(data, dataA, learnInfo, M, LForDecomp, omegaA, deltaA);

%Multiply by K_A.
MultByKNoNoiseTermA = @(x)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_re, P_r, P_c, rhoVect, newRawIndices, deltaA, ls, ms);

%Multiply by K + sigmaI.
MultByWholeK = @(x) MultByKNoNoiseTerm(x) + MultByKNoNoiseTermA(x) + sigma*x;

%Multiply by entire (K+sigmaI)^-1.
multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, CG_ITER_LIMIT);

%Now let's test matrix multiplication.
multMatByKInv = @(X) RunCGOnMatrixInitGuesser(MultByWholeK, X, PreConInv, CGErrorTol, CG_ITER_LIMIT);



visualize_phis_CG(sysInfo,obsInfo,learnInfo,'E', multMatByKInv);
visualize_phis_CG(sysInfo,obsInfo,learnInfo,'A', multMatByKInv);
toc


% tic
% learnInfo.option = 'subset';  % doesn't work for ODS
% learnInfo.Nsub = 8;
% learnInfo.sub = randsample(1:learnInfo.N,learnInfo.Nsub);
% [fval, dfval,~] = GlikSTEPreConSTESUBSET(learnInfo, learnInfo.hyp, mFORGLIK, lFORGLIK, CGErrorTol, HVAL, M, rangeOfI)
% toc
% 

return



%Visualize phis.
visualize_phis_CG(sysInfo,obsInfo,learnInfo,'E', multMatByKInv);
visualize_phis_CG(sysInfo,obsInfo,learnInfo,'A', multMatByKInv);


end