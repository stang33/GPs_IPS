function [fval, dfval,learnInfo] = GlikSTEPreConSTE(learnInfo, hyp, m, l, CGErrorTol, DerivHVal, M, rangeOfI)

%TODO: Only FM "works" rn.


STRICT_CG_ITER_LIMIT = 1000;
StrictCGErrorTol = 10^(-10);

X = learnInfo.X;
Y = learnInfo.Y;
dN = learnInfo.d*learnInfo.N*learnInfo.order;
L = length(X)/dN;
LForDecomp = L / M;
d = learnInfo.d;
N = learnInfo.N;
n = N;
D = d;

%New constant settings.
% M = obsInfo.M;
% LForDecomp = obsInfo.L;
% n = learnInfo.N;
% D = learnInfo.d;



dN1 = d*N;

%Assign the hyperparams properly.
deltaE = exp(hyp(1));
omegaE = exp(hyp(2));
deltaA = exp(hyp(3));
omegaA = exp(hyp(4));
sigma = exp(hyp(5))^2;

%If sigma is NaN, there is no noise. Use jitter factor.
if isnan(sigma)
    sigma = 10^(-6);
end


%Make preconditioner for full K matrix.
jitter = 10^(-6);
[logDetPreCon, PreConInvRaw, ~] = ConstructNystPreconNoK(learnInfo, LForDecomp, M, rangeOfI, jitter);
PreConInv = @(x) PreConInvRaw*x;


data = learnInfo.xpath_train(1:D*n,:,:);
%dataA = learnInfo.dxpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train((D*n+1):2*D*n,:,:);




%First decomp.
[U_reE, P_rE, P_cE, rhoVectE, ~, newRawIndicesE, dsE, ms, ls] = SparseDecomp(data, data, learnInfo, M, LForDecomp, omegaE, deltaE);

%Multiply by K_E.
MultByKNoNoiseTerm = @(x)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_reE, P_rE, P_cE, rhoVectE, newRawIndicesE, deltaE, ls, ms);

%Decomp for K_A.
[U_reA, P_rA, P_cA, rhoVectA, ~, newRawIndicesA, dsA, ms, ls] = SparseDecomp(data, dataA, learnInfo, M, LForDecomp, omegaA, deltaA);

%Multiply by K_A.
MultByKNoNoiseTermA = @(x)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_reA, P_rA, P_cA, rhoVectA, newRawIndicesA, deltaA, ls, ms);

%Multiply by K + sigmaI.
MultByWholeK = @(x) MultByKNoNoiseTerm(x) + MultByKNoNoiseTermA(x) + sigma*x;

%Multiply by entire (K+sigmaI)^-1.
multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, l);




%Keeping these the same from the standard Glik.
Ym = 0*Y;
Yda = 0*Y;
Ydb = 0*Y;
loga = hyp(6);
logb = hyp(7);

for i = 1:L
    for j = 1:N
        Ym((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) - (exp(loga) - exp(logb)*sum(X((dN*(i-1)+1+dN1+d*(j-1)):(dN*(i-1)+dN1+d*j)).^2))*X((dN*(i-1)+1+dN1+d*(j-1)):(dN*(i-1)+dN1+d*j));
        Yda((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) = - X((dN*(i-1)+1+dN1+d*(j-1)):(dN*(i-1)+dN1+d*j))*exp(loga);
        Ydb((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) = sum(X((dN*(i-1)+1+dN1+d*(j-1)):(dN*(i-1)+dN1+d*j)).^2)*X((dN*(i-1)+1+dN1+d*(j-1)):(dN*(i-1)+dN1+d*j))*exp(logb);

    end
end

learnInfo.Ym = Ym;

u = StandardPCG(MultByWholeK, Ym, PreConInv, StrictCGErrorTol, STRICT_CG_ITER_LIMIT);





%Store gammas for fval calculation.
gammas = zeros(l,1);

%Run standard Hutchinson with f(K) = log(K).
for i = 1 : l

    %Make a random vector of ones and zeros.
    tildez_i = 2 * randi(2,n*D*M*LForDecomp,1) - 3 * ones(n*D*M*LForDecomp,1);
    z_i = tildez_i / norm(tildez_i);

    %Do PCG and track coefs.
    [~,T] = StandardPCGWithT(MultByWholeK, z_i, PreConInv, CGErrorTol, l+1, l);

    %This is decreasing order. Doesn't matter n=8 - TODO investigate further.
    [W,lambdas] = eig(T);
    gamVal = 0;

    for j = 1 : m
        gamVal = gamVal + (W(1,j))^2 * log(lambdas(j,j));
    end

    gammas(i) = gamVal; 

end

%Some gammas might be NaN if converges exactly before l.
%Right now, we just throw out these values.
indexesNAN = isnan(gammas);
gammas(indexesNAN) = 0;

%Now we adjust l so we don't divide by extra weights we didn't use.
trueL = l - sum(indexesNAN);

%Make tau.
tau_star = logDetPreCon / 2 + (n*M*D*LForDecomp / (2*trueL)) * sum(gammas);

%Here's a good approx for fval.
fval = 1/2*Ym'*u + tau_star + log(2*pi)*(dN1*L)/2;







%Now it's dfval time.
dfval = 0*hyp;

%These are the same as Glik. No K needed.
dfval(6) = u'* Yda;
dfval(7) = u'* Ydb;



%Estimate the noise derivative.
%tau_star = JustCalculateExactTrace(multByKInv, n*D*M*LForDecomp);
tau_star = HutchinsonEst(multByKInv, l, n*D*M*LForDecomp);

%Use the trace commutative identity tr(xx') = x'x.
dfval(5) = -u' * u * sigma + tau_star * sigma;




%Now, attempt the delta deriv. No noise in K for this matrix building.
%TODO: Why not 2/delta? Is kphimatern wrong in the main code?
%MultBydKdeltaE = @(x) (2) * multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_re, P_r, P_c, rhoVect, newRawIndices, delta, ls, ms);
MultBydKdeltaE = @(x) (2) * MultByKNoNoiseTerm(x);

%Now do the trace estimation for the second term of this derivative.
MultByKInvTimesdKdeltaE = @(x) multByKInv(MultBydKdeltaE(x));

%traceEstKInvTimesdKdelta = JustCalculateExactTrace(MultByKInvTimesdKdelta, n*D*M*LForDecomp);

traceEstKInvTimesdKdeltaE = HutchinsonEst(MultByKInvTimesdKdeltaE, l, n*D*M*LForDecomp);

dfval(1) = -1*(1/2 * u' * MultBydKdeltaE(u) - 1/2 * traceEstKInvTimesdKdeltaE);



%Same thing but for A.
MultBydKdeltaA = @(x) (2) * MultByKNoNoiseTermA(x);

%Now do the trace estimation for the second term of this derivative.
MultByKInvTimesdKdeltaA = @(x) multByKInv(MultBydKdeltaA(x));

%traceEstKInvTimesdKdelta = JustCalculateExactTrace(MultByKInvTimesdKdelta, n*D*M*LForDecomp);

traceEstKInvTimesdKdeltaA = HutchinsonEst(MultByKInvTimesdKdeltaA, l, n*D*M*LForDecomp);

dfval(3) = -1*(1/2 * u' * MultBydKdeltaA(u) - 1/2 * traceEstKInvTimesdKdeltaA);





%Now for the interesting one. dKoE.
%Now let's make new rhoVects and uh, multiply by those K's.
hVal = DerivHVal;

%Make the rhoVects needed for K mult.
omegaH = omegaE + hVal;
rhoVectH = GetRhoVect(n, M, LForDecomp, omegaH, deltaE, dsE);

omegaHm = omegaE - hVal;
rhoVectHm = GetRhoVect(n, M, LForDecomp, omegaHm, deltaE, dsE);

omegaH2 = omegaE + hVal / 2;
rhoVectH2 = GetRhoVect(n, M, LForDecomp, omegaH2, deltaE, dsE);

omegaH2m = omegaE - hVal / 2;
rhoVectH2m = GetRhoVect(n, M, LForDecomp, omegaH2m, deltaE, dsE);

%Methods!
MultByKGivenRhoE = @(x,rho)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_reE, P_rE, P_cE, rho, newRawIndicesE, deltaE, ls, ms);
MultByKPartialOmegaE = @(x) omegaE * ((4/(3*hVal)) * (MultByKGivenRhoE(x,rhoVectH2) - MultByKGivenRhoE(x,rhoVectH2m)) - (1/(6*hVal)) * (MultByKGivenRhoE(x,rhoVectH) - MultByKGivenRhoE(x,rhoVectHm)));

MultByKInvTimesdKomegaE = @(x) multByKInv(MultByKPartialOmegaE(x));
%traceEstKInvTimesdKomega = HutchinsonEst(MultByKInvTimesdKomega, l, n*D*M*LForDecomp);
%traceEstKInvTimesdKomega = JustCalculateExactTrace(MultByKInvTimesdKomega, n*D*M*LForDecomp);

traceEstKInvTimesdKomegaE = HutchinsonEst(MultByKInvTimesdKomegaE, l, n*D*M*LForDecomp);

dfval(2) = -1*(1/2 * u' * MultByKPartialOmegaE(u) - 1/2 * traceEstKInvTimesdKomegaE);













%Now for the interesting one. dKoA.
%Now let's make new rhoVects and uh, multiply by those K's.
hVal = DerivHVal;

%Make the rhoVects needed for K mult.
omegaH = omegaA + hVal;
rhoVectH = GetRhoVect(n, M, LForDecomp, omegaH, deltaA, dsA);

omegaHm = omegaA - hVal;
rhoVectHm = GetRhoVect(n, M, LForDecomp, omegaHm, deltaA, dsA);

omegaH2 = omegaA + hVal / 2;
rhoVectH2 = GetRhoVect(n, M, LForDecomp, omegaH2, deltaA, dsA);

omegaH2m = omegaA - hVal / 2;
rhoVectH2m = GetRhoVect(n, M, LForDecomp, omegaH2m, deltaA, dsA);

%Methods!
MultByKGivenRhoA = @(x,rho)multByKernelGeneralNoNoiseTerm(x, sigma, n, D, M, LForDecomp, U_reA, P_rA, P_cA, rho, newRawIndicesA, deltaA, ls, ms);
MultByKPartialOmegaA = @(x) omegaA * ((4/(3*hVal)) * (MultByKGivenRhoA(x,rhoVectH2) - MultByKGivenRhoA(x,rhoVectH2m)) - (1/(6*hVal)) * (MultByKGivenRhoA(x,rhoVectH) - MultByKGivenRhoA(x,rhoVectHm)));

MultByKInvTimesdKomegaA = @(x) multByKInv(MultByKPartialOmegaA(x));
%traceEstKInvTimesdKomega = HutchinsonEst(MultByKInvTimesdKomega, l, n*D*M*LForDecomp);
%traceEstKInvTimesdKomega = JustCalculateExactTrace(MultByKInvTimesdKomega, n*D*M*LForDecomp);

traceEstKInvTimesdKomegaA = HutchinsonEst(MultByKInvTimesdKomegaA, l, n*D*M*LForDecomp);

dfval(4) = -1*(1/2 * u' * MultByKPartialOmegaA(u) - 1/2 * traceEstKInvTimesdKomegaA);




end





