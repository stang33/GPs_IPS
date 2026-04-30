    [w_SBL, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Bernoull', G_0, b_0);


    tic;
    sigma2 = var(b_0)/1e3 ;
    delta_La = 1e-10;
    lambda_init = [];
    sparsity = [];
    [w_FL,used_ids,loglik,sigma2,errbars,basis,selected,deleted,alpha,lambda] = FastLaplace(G_0,b_0, sigma2, delta_La,lambda_init,sparsity);
    time_SBL = toc;