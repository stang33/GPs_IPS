function [phi_mean,phi_cov] = phi_pos(r,learnInfo, kernel_type)
% Input: r:  1 x 1 
%        dX:  dN x 1  position data
%        kernel_type: kernel for 'E' or 'A'
% Output: K_r: 1 x dN covariance function

% (c) XXXX

X = learnInfo.X;
K = learnInfo.K;
hyp = learnInfo.hyp;
CoefM = learnInfo.CoefM;


if learnInfo.order == 1 || kernel_type == 'E'
    sigma = exp(hyp(1));
    omega = exp(hyp(2));
else
    sigma = exp(hyp(3));
    omega = exp(hyp(4));
end
    

Z=K_r(r,X,learnInfo,kernel_type);

% Ym = learnInfo.Ym;

% phi_mean = Z * pinv(K)*Ym;
phi_mean = Z * CoefM;

phi_cov = diag(cov_Matern(r,r,sigma,omega,learnInfo.v))-Z*pinv(K)*Z';

     
end





    






