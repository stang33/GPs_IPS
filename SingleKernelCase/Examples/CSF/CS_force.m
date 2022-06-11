function f = CS_force(v,d)

% input: v is the velocity data in R^dN
% F_i = a*v_i*(1-|v_i|^b)
% (c) XXXX


a = 1;  %sigma
b = 2;  %p

n =size(v,1)/d;
f=v;

for i=1:n
f((i-1)*d+1:i*d) = a*(1 - norm(v((i-1)*d+1:i*d),2)^b)*v((i-1)*d+1:i*d);
end