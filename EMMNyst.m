function [logDetPreCon, PreConInvRaw] = EMMNyst(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)

    N = learnInfo.N;
    D = learnInfo.d;

    if abs(sigma) < 10^(-6)
        sigma = jitter;
    end

    K = KE + KA + jitter * eye(M*LForDecomp*N*D);

    rank = rangeOfI;

    %Instead of randperm, we just create a custom perm to fit the EMM
    %ordering maximally well.
    ds = learnInfo.ds;
    [ind, ~] = FPS(ds, rangeOfI); %i pts from ds

    perm = zeros(1,N*D*M*LForDecomp);

    %We need to extract the right particles from the pairs.
    %For each index:

    pairs = learnInfo.pairs;
    fullInd = zeros(1,2*rangeOfI);

    for i = 1 : rangeOfI

        curDS = ind(i);

        curM = ceil(curDS / ((N*(N-1))/2*LForDecomp));
        curL = ceil(curDS / ((N*(N-1))/2)) - (curM-1) * LForDecomp;

        %This is the place in the ordered permutation list.
        place = mod(curDS - 1, ((N*(N-1))/2)  ) + 1;
        curX = pairs(1,1,place);
        curY = pairs(1,2,place);
        
        fullInd(1,2*i-1) = (curM-1) * (N*D*LForDecomp) + (curL-1) * (N*D) + curX;
        fullInd(1,2*i) = (curM-1) * (N*D*LForDecomp) + (curL-1) * (N*D) + curY;

    end
    

    perm(1:2*rangeOfI) = fullInd;

    curIndex = 2*rangeOfI + 1;
    for i = 1 : N*D*M*LForDecomp
        if ~ismember(i,perm)
            perm(1,curIndex) = i;
            curIndex = curIndex + 1;
        end
    end

    P = eye(N*D*M*LForDecomp);
    P = P(perm,:);


    %Run Nyst code.

    %Initialize
    Aplus = zeros(rank, rank);
    
    %Permute rows and columns to get the columns we want to use.
    PTKP = K(perm, perm);
    A = PTKP(1:rank, 1:rank);
    
    [U,S,~] = svd(A,"vector");
    
    for i = 1 : rank
        
        singval = S(i);
        col = U(:,i);
        Aplus = Aplus + (1.0 / singval) * (col * col.');
    
    end
    
    R = PTKP(1:rank,:);

    %Get output.
    lognoise = learnInfo.hyp(5); %noise of Y

    %Woodbury inversion.
    PreConInvRaw = (1 / exp(lognoise)^2)*eye(N*D*M*LForDecomp) - (1 / exp(lognoise)^4)*P'*R'*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R * R')))*R*P;

    %Symmetrize.
    mat = pinv(Aplus) + 1/sigma * (R * R');
    symmat = (mat + mat') / 2;

    logDetPreCon = trace(logm(symmat)) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma);
    

end




