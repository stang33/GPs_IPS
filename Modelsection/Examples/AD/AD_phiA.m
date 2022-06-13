function phi = AD_phiA(r)
% function phi = AD_phiA(r, tau, U_prime)

% (C)
tau = 0.1;

phi = tau * (1+r.^2).^(-0.5);
end