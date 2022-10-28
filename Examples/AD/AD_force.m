function f = AD_force(v,d)

% input: v is the velocity data in R^dN
% F_i = ...

n =size(v,1)/d;
f=v;

for i=1:n
f((i-1)*d+1:i*d) = 0*v((i-1)*d+1:i*d);
end