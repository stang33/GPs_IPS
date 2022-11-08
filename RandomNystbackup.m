function [logDetPreCon, PreConInvRaw] = RandomNystbackup(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)

    N = learnInfo.N;
    D = learnInfo.d;

    if abs(sigma) < 10^(-6)
        sigma = jitter;
    end

    K = KE + KA + jitter * eye(M*LForDecomp*N*D);

    K = (K + K') / 2;

    rank = rangeOfI;

    perm = randperm(N*D*M*LForDecomp);
    P = eye(N*D*M*LForDecomp);
    P = P(perm,:);


    %Initialize
    Aplus = zeros(rank, rank);
    
    %Permute rows and columns to get the columns we want to use.
    PTKP = K(perm, perm);

    A = PTKP(1:rank, 1:rank);
    
    [U,S,~] = svd(A,"vector");
    
    for i = 1 : rank
        
        singval = S(i);

        if singval > 0
            col = U(:,i);
            Aplus = Aplus + (1.0 / singval) * (col * col.');
        end

    end

    
    R = PTKP(1:rank,:);

    lognoise = learnInfo.hyp(5); %noise of Y
    PreConInvRaw = (1 / exp(lognoise)^2)*eye(N*D*M*LForDecomp) - (1 / exp(lognoise)^4)*P'*R'*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R * R')))*R*P;
  
    mat = pinv(Aplus) + 1/sigma * (R * R');
    symmat = (mat + mat') / 2;
    
    logDetPreCon = trace(logm(symmat)) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma);
        
    
    end



