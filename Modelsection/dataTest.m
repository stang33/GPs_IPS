
addpaths;
load("FM10L3M3");
%load("AD10L3M4.mat");
%load("FM10nu52");

M = obsInfo.M;
LForDecomp = obsInfo.L;
n = learnInfo.N;
D = learnInfo.d;




return



%  load("CSF10L3M3.mat");
%  obsInfo.L = 3;

%learnInfo.hyp


rangeOfI = 6;
jitter = 10^(-6);

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


%Adjustable parameters.
CGErrorTol = 10^(-5);
CG_ITER_LIMIT = 40;

%Decomp for K_E.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train(D*n+1:2*D*n,:,:);
% dataAA = learnInfo.dxpath_train(1:D*n,:,:);
% dataAAA = learnInfo.dxpath_train(D*n+1:2*D*n,:,:);





trueK = learnInfo.K;

[U_s, R_s, UnsortU, sortingInd] = TotalDecompForDebug(data, omegaE, deltaE, n, D, M, LForDecomp);
U_s*R_s*U_s' - trueK
norm(U_s*R_s*U_s' - trueK)


return





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

%Compute the Ym product once.
learnInfo.invKTimesYm = multByKInv(learnInfo.Ym);














return
%TODO: Add Preconditioner.
%PreConInv = @(x) x;
[logDetPreCon, PreConInvMat, PreCon] = ConstructNystPreconNoK(learnInfo, LForDecomp, M, rangeOfI, jitter);
PreConInv = @(x) PreConInvMat * x;




trueK = learnInfo.K;
cond(trueK)
cond(PreConInv(trueK))
norm(trueK - PreCon)
norm(pinv(trueK) - PreConInvMat)
norm(pinv(trueK) - pinv(PreCon))
norm(PreConInvMat - pinv(PreCon))

%pinv(trueK) - pinv(PreCon)
% mat = PreConInvMat - pinv(PreCon)
% mat(1,:)
%PreConInvMat - pinv(PreCon)
return



%Visualize phis.
visualize_phis_CG(sysInfo,obsInfo,learnInfo,'E', multMatByKInv);
visualize_phis_CG(sysInfo,obsInfo,learnInfo,'A', multMatByKInv);


end