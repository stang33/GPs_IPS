function y = cov_Matern(x_1,x_2,sigma,omega,v)

% Matern covariance function with nu = d/2 and isotropic distance measure. For
% d=1 the function is also known as the exponential covariance function or the 
% Ornstein-Uhlenbeck covariance in 1d. The covariance function is:
%
%   k(x,z) = f( sqrt(d)*r ) * exp(-sqrt(d)*r)
%
% with f(t)=1 for d=1, f(t)=1+t for d=3, f(t)=1+t+t^2/3 for d=5 and
%      f(t)=1+t+2*t^2/5+t^3/15  for d=7.
%
% The covariance function can also be expressed for non-integer d.
%
%   k(x,z) = 2^(1-nu)/gamma(nu) * s^nu * besselk(nu,s), s = sqrt(2*nu)*r
%
% Note that for d->oo the covariance converges to the squared exponential.
%
% Here r is the Mahalanobis distance sqrt(maha(x,z)). The function takes a
% "mode" parameter, which specifies precisely the Mahalanobis distance used, see
% covMaha. The function returns either the number of hyperparameters (with less
% than 3 input arguments) or it returns a covariance matrix and (optionally) a
% derivative function.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-05-29.
%
% See also cov/covSE.m, cov/covMaha.m.

r = abs(x_1-x_2).*omega;

switch v
    case 1/2
        
        y = sigma^2.*exp(-r);
    
    case 3/2
        
        y = sigma^2*(1+sqrt(3).*r).*exp(-sqrt(3).*r);
        
    case 5/2
        y = sigma^2*(1+sqrt(5).*r+5.*r.^2/3).*exp(-sqrt(5).*r);
        
        
     case 7/2
         
        y = sigma^2*(1+sqrt(7).*r+0.4*r.^2+ r.^3/15).*exp(-sqrt(7).*r);
        
        
        
end

end

