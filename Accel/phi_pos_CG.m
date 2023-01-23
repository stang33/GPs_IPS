function [phi_mean,phi_cov] = phi_pos_CG(r,learnInfo, kernel_type, multKInvByMatrix)
% Input: r:  1 x 1 
%        dX:  dN x 1  position data
%        kernel_type: kernel for 'E' or 'A'
% Output: K_r: 1 x dN covariance function

% (c) XXXX

X = learnInfo.X;
hyp = learnInfo.hyp;

if learnInfo.order == 1 || kernel_type == 'E'
    sigma = exp(hyp(1));
    omega = exp(hyp(2));
else
    sigma = exp(hyp(3));
    omega = exp(hyp(4));
end
    

Z=K_r(r,X,learnInfo,kernel_type);

phi_mean = Z * learnInfo.invKTimesYm;

phi_cov = diag(cov_Matern(r,r,sigma,omega,learnInfo.v))-Z*multKInvByMatrix(Z');

     
end





    






