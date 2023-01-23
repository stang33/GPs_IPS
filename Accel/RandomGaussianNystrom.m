function [logDetPreCon, PreConInvRaw] = RandomGaussianNystrom(rank, jitter, KE, KA, sigma)

    %Randomized Gaussian Nystrom approximation.

    if abs(sigma) < 10^(-6)
        sigma = jitter;
    end
    
    %Construct K.
    K = KE + KA;
    K = (K + K') / 2;

    %Follow algorithm.
    Omega = randn(size(K,1),rank);
    Omega = qr(Omega);
    Y= K*Omega;
    v = eps(norm(Y,'fro'));
    Yv = Y+v*Omega;
    C = chol(Omega'*Yv);
    B = Yv/C;
    [U,Sigma,~] = svd(B,0);
    Lambda = max(0,Sigma^2-v*eye(size(Sigma^2,1)));

    %Return the raw inverse matrix.
    PreConInvRaw = (Lambda(end)+sigma)*U*pinv(Lambda+sigma*eye(size(U,2)))*U'+eye(size(K,1))-U*U';
    v= reshape(Lambda(Lambda>1e-6),[],1);  

    %Return the log det.
    logDetPreCon = sum(sum(log(1./(v+sigma))));
    
end



