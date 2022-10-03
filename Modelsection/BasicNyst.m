function [logDetPreCon, PreConInvRaw] = BasicNyst(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)

    N = learnInfo.N;
    D = learnInfo.d;
    %Let's say rangeOfI is the number of rows in K?
    %PreConInvRaw = eye(M*LForDecomp*N*D);

    if abs(sigma) < 10^(-8)
        sigma = jitter;
    end

    %K = KE + KA + sigma * eye(M*LForDecomp*N*D);
    K = KE + KA + jitter * eye(M*LForDecomp*N*D);

    rank = rangeOfI;


    %Run Nyst code.

    %Initialize
    Aplus = zeros(rank, rank);
    
    %Permute rows and columns to get the columns we want to use.
    %PTKP = K(1:rank, 1:rank);
    
    A = K(1:rank, 1:rank);
    
    [U,S,~] = svd(A,"vector");
    
    for i = 1 : rank
        
        singval = S(i);

        %Try a cutoff?
        %if singval > 10^(-4)
            col = U(:,i);
            Aplus = Aplus + (1.0 / singval) * (col * col.');
        %end
    
    end
    
    R = K(:,1:rank);
    
%Diff = K - R*Aplus*R'


    %Get output.
    lognoise = learnInfo.hyp(5); %noise of Y
    PreConInvRaw = (1 / exp(lognoise)^2)*eye(M*LForDecomp*N*D) - (1 / exp(lognoise)^4)*R*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R' * R)))*R';
    
    
    logDetPreCon = trace(logm(pinv(Aplus) + 1/sigma * (R' * R))) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma);
    

% trace(logm(pinv(Aplus) + 1/sigma * (R' * R)))
% trace(logm(Aplus))
% D*M*N*LForDecomp * log(sigma)
% 
% hi = 11111111


% cond(KE + KA + sigma * eye(M*LForDecomp*N*D))
% cond(PreConInvRaw*(KE + KA + sigma * eye(M*LForDecomp*N*D)))



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


