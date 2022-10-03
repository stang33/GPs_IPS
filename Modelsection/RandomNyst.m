function [logDetPreCon, PreConInvRaw] = RandomNyst(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)

    N = learnInfo.N;
    D = learnInfo.d;
    %Let's say rangeOfI is the number of rows in K?
    %PreConInvRaw = eye(M*LForDecomp*N*D);

    if abs(sigma) < 10^(-6)
        sigma = jitter;
    end

    K = KE + KA + jitter * eye(M*LForDecomp*N*D);

    rank = rangeOfI;



    %zerop = 1:N*D*M*LForDecomp;
    perm = randperm(N*D*M*LForDecomp);
    %invperm(perm) = zerop(1:N*D*M*LForDecomp);
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
    

%     P = P(perm,:);
%     K = K;
%     P2 = eye(N*D*M*LForDecomp);
%     P2 = P2(invperm,:);
%     NystK = (R' * Aplus * R) *  P2';
%     NystK2 = R' * Aplus * R;
%     Diff = K - P'*NystK2*P;
%     Diff(1:20,1:20)
%     Diff2 = PTKP - NystK2;
% Diff3 = K - NystK2(invperm,invperm);
% Diff3(1:10,1:10)
% Diff2(1:10,1:10)

    %Get output.
    
    %logDetPreCon = trace(logm(pinv(Aplus) + 1/sigma * (R' * R))) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma);


    lognoise = learnInfo.hyp(5); %noise of Y
    %nonPerm_invKPlusSigma2 = (1 / exp(lognoise)^2)*P' - (1 / exp(lognoise)^4)*P'*R'*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R * P' * R')))*R*P';
    PreConInvRaw = (1 / exp(lognoise)^2)*eye(N*D*M*LForDecomp) - (1 / exp(lognoise)^4)*P'*R'*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R * R')))*R*P;
    %PreConInvRaw = nonPerm_invKPlusSigma2(invperm, invperm);

% no = cond(K+ sigma * eye(M*LForDecomp*N*D))
% OUT2 = cond(PreConInvRaw * (K + sigma * eye(M*LForDecomp*N*D)))
% 
% TRUEDET = trace(logm(pinv(PreConInvRaw)))
% TRUEDET2 = log(det(pinv(PreConInvRaw)))
% APPDET = trace(logm(pinv(Aplus) + 1/sigma * (R * R'))) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma)
logDetPreCon = trace(logm(pinv(Aplus) + 1/sigma * (R * R'))) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma);
    

end





% %Finally, compare with random.
% 
% p = randperm(90);
% [~,invperm] = sort(p);
% 
% TrueKJitter = TrueK + 10^(-6)*eye(90);
% K = TrueKJitter(p,p);
% 
% 
% 
% 
% 
% rank = pts;
% Aplus = zeros(rank, rank);
% %jitter added prior
% A = K(1:rank, 1:rank);% + 10^(-6)*eye(rank,rank);
% 
% 
% [U,S,~] = svd(A,"vector");
% 
% for i = 1 : rank
%     
%     singval = S(i);
%     col = U(:,i);
%     Aplus = Aplus + (1.0 / singval) * (col * col.');
% 
% end
% 
% R = K(1:rank,:);
% 
% P = eye(90);
% P = P(p,p);
% 
% 
% %Get the inverse.
% lognoise = learnInfo.hyp(5); %noise of Y
% nonPerm_invKPlusSigma2 = (1 / exp(lognoise)^2)*P' - (1 / exp(lognoise)^4)*P'*R'*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R * P' * R')))*R*P';
% PreCInv3 = nonPerm_invKPlusSigma2(invperm, invperm);
% 
% 
% 
% 
% TrueKNoise = TrueK + sigma*eye(90);
% 
% cond(TrueKNoise)
% cond(PreCInv3 * TrueKNoise)





% 
% function [R, Aplus] = Nyst(K, rank, perm)
% 
%     %Initialize
%     Aplus = zeros(rank, rank);
%     
%     %Permute rows and columns to get the columns we want to use.
%     PTKP = K(perm, perm);
%     
%     A = PTKP(1:rank, 1:rank);
%     
%     [U,S,~] = svd(A,"vector");
%     
%     for i = 1 : rank
%         
%         singval = S(i);
%         col = U(:,i);
%         Aplus = Aplus + (1.0 / singval) * (col * col.');
%     
%     end
%     
%     R = PTKP(:,1:rank);
%     
%     %R = R(:,perm);
% 
% end
% 
% 
% 
% 
% 
% function output = ConstructNystPrecon(learnInfo, rank)
% 
%     n = learnInfo.N;
%     D = learnInfo.d;
% 
%     zerop = 1:n*D;
%     perm = randperm(n*D);
%     invperm(perm) = zerop(1:n*D);
%     P = eye(n*D);
%     P = P(perm,perm);
%     
%     
%     [R, Aplus] = Nyst(learnInfo.K, rank, perm);
%     
%     
%     %Get the inverse. TODO Make method of this for fast CG.
%     lognoise = learnInfo.hyp(5); %noise of Y
%     nonPerm_invKPlusSigma2 = (1 / exp(lognoise)^2)*P' - (1 / exp(lognoise)^4)*P'*R*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R' * P' * R)))*R'*P';
%     PreCInv = nonPerm_invKPlusSigma2(invperm, invperm);
%     output = PreCInv;
% 
% 
% end
% 
%         


