function [psi, dpsi] = mcos(r, omega, p)

psi = cos(r*omega+p);
dpsi = -omega*sin(r*omega+p);

return