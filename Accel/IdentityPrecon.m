function [logDetPreCon, PreConInvRaw] = IdentityPrecon(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)

    N = learnInfo.N;
    D = learnInfo.d;
    logDetPreCon = 0;
    PreConInvRaw = eye(M*LForDecomp*N*D);