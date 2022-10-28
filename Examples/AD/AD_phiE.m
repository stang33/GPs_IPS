function phi = AD_phiE(r)
% function phi = AD_phiE(r, s, tau, U_prime, U_dprime)

% (C) M. Zhong

tau = 0.1;

phi = tau * (1+r).^(-2.5) + (1+r).^(-0.5);

end