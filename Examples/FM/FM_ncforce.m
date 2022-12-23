function f = FM_ncforce(v,d)

% input: v is the velocity data in R^dN
% (c) XXXX


alpha = 1.5;
beta = 0.5;
n =size(v,1)/d;
f=v;

for i=1:n
f((i-1)*d+1:i*d) = (alpha - beta*norm(v((i-1)*d+1:i*d),2)^2)*v((i-1)*d+1:i*d);
end