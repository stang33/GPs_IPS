function [psi, dpsi] = mexp(r, omegasq, weight)

psi = weight.*exp(-r.*r*omegasq);
dpsi = -2*r*weight*omegasq.*psi;

return