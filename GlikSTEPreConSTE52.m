function [fval, dfval,learnInfo] = GlikSTEPreConSTE52(learnInfo, hyp, m, l, CGErrorTol, ~, M, rangeOfI, PreconMethod)

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
[U_sE, R_sE, ~, ~] = TotalDecompForDebug52(data, data, omegaE, deltaE, n, D, M, LForDecomp);
KE = U_sE * R_sE * U_sE';
MultByKE = @(x) KE * x;

[U_sA, R_sA, ~, ~] = TotalDecompForDebug52(data, dataA, omegaA, deltaA, n, D, M, LForDecomp);
KA = U_sA * R_sA * U_sA';
MultByKA = @(x) KA * x;

%Make preconditioner for full K matrix.
jitter = 10^(-6);
[logDetPreCon, PreConInvRaw] = PreconMethod(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma);
PreConInv = @(x) PreConInvRaw*x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Let me try to understand
% I am talking about RandomNyst.m
% Q1: Is P for K or K+sigmaI? 
% Q2: PreConInvRaw = inv(P)?
% Q3: logDetPreCon = log(det(P))?
% 
% Q1: P is an approximation of only K, not K+sigmaI. This is because
%     our goal is to construct a preconditioner that is easy to invert.
%     We hope that P if approximates K, then P+sigmaI will approximate K+sigmaI.
%     This allows us to invert P+sigmaI quickly, as the Woodbury inversion
%     formula allows us to only invert explicitly matrices of reduced rank,
%     namely the rank of our preconditioner.
%
% Q2: PreConInvRaw = inv(P + sigmaI)
%
% Q3: logDetPreCon = log(det(P + sigmaI)). By taking the log of the Woodbury
%     formula, we can avoid taking the trace log of the entire matrix, and
%     again only evaluate matrix logarithms up to the rank.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%Multiply by K + sigmaI.
MultByWholeK = @(x) MultByKE(x) + MultByKA(x) + sigma*x;

%Multiply by entire (K+sigmaI)^-1.
multByKInv = @(x) StandardPCG(MultByWholeK, x, PreConInv, CGErrorTol, l);


%These are the same from the standard Glik.
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

    %This is decreasing order.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I can not understand two things
% Q1: why it is 2*trueL? according to formulas, is it the number of vector
% used in probing vectors?
% Q2: The weights coming from the matrix T is to compute log(Det(K+sigma)
%  I think if you have  computed logDetPreCon=log(det(P)), the right thing to compute
%  should be logdet(inv(P)*(K+sigmaI)), are you sure this is
%  (n*M*D*LForDecomp / (2*trueL)) * sum(gammas);?
%
% Q1: We need 1/2 * fval. The trueL issue is an interesting one. There are
% times that we converge to numerical roundoff error in our PCG BEFORE we have
% completed the desired number of PCG runs. This will be reflected in
% Lanczos weights that are 0. Thus we will have a gamma value that is NaN.
% We disregard these, and generally we must just be careful to not request
% too high of a number of Lanczos weights that we consistently hit
% numerical zero before completiton, or we may accidentally have
% sum(gammas) = 0 and trueL = 0, which means fval is NaN. If this ever
% happens, huge errors are reported by the compiler, and I believe I have
% tuned the values well enough to consistently have very very few gamma
% values that error out to NaN.
%
% Q2: The term (n*M*D*LForDecomp / (2*trueL)) * sum(gammas) approximates
% the log det[(P + sigmaI)^(-1) (K + sigmaI)]. This formula comes from
% the Lanczos algorithm. We then add the log(det(P+sigmaI)) and rescale by
% our constants to ensure we match the form of the actual maximum
% likelihood function when using a kernel of our specified size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%Here's a good approx for fval.
fval = 1/2*Ym'*u + tau_star + log(2*pi)*(dN1*L)/2;




%Now it's dfval computation.
dfval = 0*hyp;

%These are the same as Glik. No K needed.
dfval(6) = u'* Yda;
dfval(7) = u'* Ydb;


%Estimate the noise derivative.
tau_star = HutchPlusPlus(multByKInv, l, n*D*M*LForDecomp);

%Use the trace commutative identity tr(xx') = x'x.
dfval(5) = -u' * u * sigma + tau_star * sigma;


%Now, attempt the delta deriv.
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




[U_sE, R_sE, ~, ~] = TotalDecompForDebug52Partial(data, data, omegaE, deltaE, n, D, M, LForDecomp);
KEPartial = U_sE * R_sE * U_sE';
MultKEPartial = @(x) KEPartial * x;

MultByKInvTimesdKomegaE = @(x) multByKInv(MultKEPartial(x));
traceEstKInvTimesdKomegaE = HutchPlusPlus(MultByKInvTimesdKomegaE, l, n*D*M*LForDecomp);
dfval(2) = -1*(1/2 * u' * MultKEPartial(u) - 1/2 * traceEstKInvTimesdKomegaE);



[U_sA, R_sA, ~, ~] = TotalDecompForDebug52Partial(data, dataA, omegaA, deltaA, n, D, M, LForDecomp);
KAPartial = U_sA * R_sA * U_sA';
MultKAPartial = @(x) KAPartial * x;

MultByKInvTimesdKomegaA = @(x) multByKInv(MultKAPartial(x));
traceEstKInvTimesdKomegaA = HutchPlusPlus(MultByKInvTimesdKomegaA, l, n*D*M*LForDecomp);
dfval(4) = -1*(1/2 * u' * MultKAPartial(u) - 1/2 * traceEstKInvTimesdKomegaA);




end





