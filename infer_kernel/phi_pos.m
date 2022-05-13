function [phi_mean,phi_cov] = phi_pos(r,learnInfo)
% Input: r:  1 x 1 
%        dX:  dN x 1  position data
% Output: K_r: 1 x dN covariance function

% (c) XXXX

X = learnInfo.X;
invK = learnInfo.invK;
invKprodYm = learnInfo.invKprodYm;
hyp = learnInfo.hyp;
sigma = exp(hyp(1));
omega = exp(hyp(2));

Z=K_r(r,X,learnInfo);

phi_mean = Z * invKprodYm;

phi_cov = diag(cov_Matern(r,r,sigma,omega,learnInfo.v))-Z*invK*Z';

     
end





    






