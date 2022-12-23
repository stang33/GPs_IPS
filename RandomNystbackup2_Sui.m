function [logDetPreCon, PreConInvRaw] = RandomNystbackup2_Sui(learnInfo, LForDecomp, M, rangeOfI, jitter, KE, KA, sigma)
% 
%     N = learnInfo.N;
%     D = learnInfo.d;
% 
%     if abs(sigma) < 10^(-6)
%         sigma = jitter;
%     end
% 
%     %K is our full kernel matrix, K_phiE + K_phiA.
%     K = KE + KA + jitter * eye(M*LForDecomp*N*D);
% 
%     %Symmetrize.
%     K = (K + K') / 2;
% 
%     %Rank parameter. rank = s. Number of columns subsampled.
%     rank = rangeOfI;
% 
%     %Choose random permutation, make a permuation matrix.
%     perm = randperm(N*D*M*LForDecomp);
%     P = eye(N*D*M*LForDecomp);
%     P = P(perm,:);
% 
% 
%     %Initialize
%     Aplus = zeros(rank, rank);
%     
%     %Permute rows and columns to get the columns we want to use.
%     PTKP = K(perm, perm);
% 
%     A = PTKP(1:rank, 1:rank);
%     
%     [U,S,~] = svd(A);
%     
%     for i = 1 : rank
%         
%         singval = S(i);
% 
%         if singval > 0
%             col = U(:,i);
%             Aplus = Aplus + (1.0 / singval) * (col * col.');
%         end
% 
%     end
% 
%     
%     R = PTKP(1:rank,:);
% 
%     lognoise = learnInfo.hyp(5); %noise of Y
%     PreConInvRaw = (1 / exp(lognoise)^2)*eye(N*D*M*LForDecomp) - (1 / exp(lognoise)^4)*P'*R'*pinv((pinv(Aplus) + (1/exp(lognoise)^2) * (R * R')))*R*P;
%   
% %     mat = pinv(Aplus) + 1/sigma * (R * R');
% %     symmat = (mat + mat') / 2;
%     
%      mat1 = pinv(Aplus);
%     symmat = mat1 + 1/sigma * (R * R');
%     %symmat = (mat + mat') / 2;
%     [~,S,~]=svd(Aplus);
%     v= reshape(S(S>1e-10),[],1);
%     
%    % logDetPreCon = trace(logm(symmat)) + trace(logm(Aplus)) + D*M*N*LForDecomp * log(sigma);
%   logDetPreCon = trace(logm(symmat)) +sum(log(v))  + D*M*N*LForDecomp * log(sigma);
%     

%randomized nystrom approximation
%Input: 
    N = learnInfo.N;
    D = learnInfo.d;

    if abs(sigma) < 10^(-6)
        sigma = jitter;
    end
    
    jitter =0;
    K = KE + KA + jitter * eye(M*LForDecomp*N*D);

    K = (K + K') / 2;

    rank = rangeOfI;

    Omega = randn(size(K,1),rank);
    Omega = qr(Omega);


    
    Y= K*Omega;
    v = eps(norm(Y,'fro'));
    Yv = Y+v*Omega;
    C = chol(Omega'*Yv);
    B = Yv/C;
    [U,Sigma,~] = svd(B,0);


    Lambda = max(0,Sigma^2-v*eye(size(Sigma^2,1)));

    %Initialize
   % Kplus = Omega'*K*Omega;
    
    %Khat = K*Omega*pinv(Kplus)*Omega'*K';% nystrom approximation of K
    
    %[U,Lambda,~] = svd(Khat,"econ");
    
    
    % construct conditioner P
    
   %sigma = learnInfo.hyp(5); %noise of Y

    
    %P = 1./(Lambda(end)+sigma)*U*(Lambda+sigma*eye(size(U,2)))*U'+eye(size(K,1))-U*U';
    % comput log (det (Khat))
 

   %save('Pmat','P','Khat');
    
    %logDetPreCon = trace(logm(P)) + trace(logm(Pinv(P))) + D*M*N*LForDecomp * log(sigma);
    PreConInvRaw = (Lambda(end)+sigma)*U*pinv(Lambda+sigma*eye(size(U,2)))*U'+eye(size(K,1))-U*U'; % inverse of P
    % PreConInvRaw = eye(size(K,1));
        
   %logDetPreCon = trace(logm(K+sigma*eye(size(K,1))));  %+ D*M*N*LForDecomp * log(sigma
    v= reshape(Lambda(Lambda>1e-6),[],1);  



    logDetPreCon = sum(sum(log(1./(v+sigma)))); % log(P+sigma
    
    save('P.mat','logDetPreCon','Lambda','v','sigma');
    end



