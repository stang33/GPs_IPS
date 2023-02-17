function [fval, dfval,learnInfo] = GlikGeneral(learnInfo, hyp, m, l, CGErrorTol, M, rank, PreconMethod, LearnBasisHyps)

%Set ultra-strict parameters for the often used u vector.
STRICT_CG_ITER_LIMIT = 1000;
StrictCGErrorTol = 10^(-10);

%Assign needed quantites from data.
X = learnInfo.X;
Y = learnInfo.Y;
dN = learnInfo.d*learnInfo.N*learnInfo.order;
L = length(X)/dN;
LForDecomp = L / M;
d = learnInfo.d;
N = learnInfo.N;
n = N;
D = d;
dN1 = d*N;
Nsub = N;
sub = 1:learnInfo.N;
dNs = d*Nsub;
name = learnInfo.name;

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

%Reference data.
data = learnInfo.xpath_train(1:D*n,:,:);
dataA = learnInfo.xpath_train((D*n+1):2*D*n,:,:);


%Create the kernel.
KE = ConstructKernelMatrix(data, data, omegaE, deltaE, n, D, M, LForDecomp, learnInfo.v, false);
MultByKE = @(x) KE * x;

KA = ConstructKernelMatrix(data, dataA, omegaA, deltaA, n, D, M, LForDecomp, learnInfo.v, false);
MultByKA = @(x) KA * x;

%Make preconditioner for full K matrix.
jitter = 10^(-6);
[logDetPreCon, PreConInvRaw] = PreconMethod(rank, jitter, KE, KA, sigma);
PreConInv = @(x) PreConInvRaw*x;

%Multiply by K + sigmaI.
MultByWholeK = @(x) MultByKE(x) + MultByKA(x) + sigma*x;

%Multiply by entire (K+sigmaI)^-1.
multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, l, 0);

Ym = zeros(dNs*L,1);
Yda = zeros(dNs*L,1);
Ydb = zeros(dNs*L,1);

if strcmp("CSF", name)

    loga = hyp(6); % parameter in noncollective force
    logb = hyp(7);
    for i = 1:L
        for j = 1:Nsub
            Ym((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j))) - (exp(loga) * X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*(1- norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb)));
            Yda((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = - X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*(1- norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb))*exp(loga);
            Ydb((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = exp(loga) * X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb)*log(norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2))*exp(logb);
        end
    end
    learnInfo.Ym = Ym;

elseif strcmp(name,"AD")

    for i = 1:L
        for j = 1:Nsub
            Ym((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j)));
        end
    end
    learnInfo.Ym = Ym;

elseif strcmp(name,"FM")

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

end



u = StandardPCG(MultByWholeK, Ym, PreConInv, StrictCGErrorTol, STRICT_CG_ITER_LIMIT, 0);


%Run Stochastic Lanczos.
%Store gammas for fval calculation.
gammas = zeros(l,1);

%Run standard Hutchinson with f(K) = log(K).
for i = 1 : l

    %Make a random vector of ones and zeros.
    tildez_i = 2 * randi(2,n*D*M*LForDecomp,1) - 3 * ones(n*D*M*LForDecomp,1);
    z_i = tildez_i / norm(tildez_i);

    %Do PCG and track coefs.
    [~,T] = StandardPCGWithT(MultByWholeK, z_i, PreConInv, CGErrorTol, m+1, m);

    [W,lambdas] = eig(T);
    gamVal = 0;

    for j = 1 : m
        if abs(lambdas(j,j)) > 10^(-8) 
            gamVal = gamVal + (W(1,j))^2 * log(lambdas(j,j));
        end
    end

    gammas(i) = gamVal; 

end

%Some gammas might be NaN if converges exactly before l.
%We just throw out these values.
indexesNAN = isnan(gammas);
gammas(indexesNAN) = 0;

%Now we adjust l so we don't divide by extra weights we didn't use.
trueL = l - sum(indexesNAN);

%Make tau.
tau_star = logDetPreCon / 2 + (n*M*D*LForDecomp / (2*trueL)) * sum(gammas);

%Computed approx for fval.
fval = 1/2*Ym'*u + tau_star + log(2*pi)*(dN1*L)/2;




%Now, dfval computation.
dfval = 0*hyp;

if strcmp(name,"FM") || strcmp(name,"CSF")
    dfval(6) = u'* Yda;
    dfval(7) = u'* Ydb;
end

%Estimate the noise derivative.
tau_star = HutchPlusPlus(multByKInv, l, n*D*M*LForDecomp);

%Use the trace commutative identity tr(xx') = x'x.
dfval(5) = -u' * u * sigma + tau_star * sigma;

if LearnBasisHyps

    MultBydKdeltaE = @(x) (2) * MultByKE(x);
    
    %Now do the trace estimation for the second term of this derivative.
    MultByKInvTimesdKdeltaE = @(x) multByKInv(MultBydKdeltaE(x));
    
    traceEstKInvTimesdKdeltaE = HutchPlusPlus(MultByKInvTimesdKdeltaE, l, n*D*M*LForDecomp);
    
    dfval(1) = -1*(1/2 * u' * MultBydKdeltaE(u) - 1/2 * traceEstKInvTimesdKdeltaE);
    
    
    
    %Same thing but for A.
    MultBydKdeltaA = @(x) (2) * MultByKA(x);
    
    %Now do the trace estimation for the second term of this derivative.
    MultByKInvTimesdKdeltaA = @(x) multByKInv(MultBydKdeltaA(x));
    
    traceEstKInvTimesdKdeltaA = HutchPlusPlus(MultByKInvTimesdKdeltaA, l, n*D*M*LForDecomp);
    
    dfval(3) = -1*(1/2 * u' * MultBydKdeltaA(u) - 1/2 * traceEstKInvTimesdKdeltaA);
    
    
    
    
    KEPartial = ConstructKernelMatrix(data, data, omegaE, deltaE, n, D, M, LForDecomp, learnInfo.v, true);
    MultKEPartial = @(x) KEPartial * x;
    
    MultByKInvTimesdKomegaE = @(x) multByKInv(MultKEPartial(x));
    traceEstKInvTimesdKomegaE = HutchPlusPlus(MultByKInvTimesdKomegaE, l, n*D*M*LForDecomp);
    dfval(2) = -1*(1/2 * u' * MultKEPartial(u) - 1/2 * traceEstKInvTimesdKomegaE);
    
    
    
    KAPartial = ConstructKernelMatrix(data, dataA, omegaA, deltaA, n, D, M, LForDecomp, learnInfo.v, true);
    MultKAPartial = @(x) KAPartial * x;
    
    MultByKInvTimesdKomegaA = @(x) multByKInv(MultKAPartial(x));
    traceEstKInvTimesdKomegaA = HutchPlusPlus(MultByKInvTimesdKomegaA, l, n*D*M*LForDecomp);
    dfval(4) = -1*(1/2 * u' * MultKAPartial(u) - 1/2 * traceEstKInvTimesdKomegaA);


end




end





