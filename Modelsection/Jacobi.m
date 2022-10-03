function [logDetPreCon, PreConInvRaw] = Jacobi(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)

    N = learnInfo.N;
    D = learnInfo.d;
    %Let's say rangeOfI is the number of rows in K?
    %PreConInvRaw = eye(M*LForDecomp*N*D);

    if abs(sigma) < 10^(-8)
        sigma = jitter;
    end

    K = KE + KA + sigma * eye(M*LForDecomp*N*D);

    PreConInvRaw = zeros(M*LForDecomp*N*D);
    logDetPreCon = 1;

    for i = 1 : N*D*M*LForDecomp

        PreConInvRaw(i,i) = 1 / K(i,i);
        logDetPreCon = logDetPreCon * K(i,i);

    end

    
    logDetPreCon = log(logDetPreCon);   
    

end



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


