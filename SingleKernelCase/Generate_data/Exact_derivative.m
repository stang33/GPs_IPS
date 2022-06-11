function [dx] =Exact_derivative(x,myode)
% compute the exact derivetives from trajectory data
%
% Input - x: trajectory data Nd x L
%
%       - myode: Nd componentwise ode equations
%
% Output 
%        -dx: the exaxt velocities at the same time instances with x 
%
% (c) XXXX

dx = x;
L  = size(x,2);
M  = size(x,3);

for j = 1:L
    for i = 1:M
        dx(:,j,i) = myode(0,x(:,j,i)); % this is autonomous ode
    end
end

end