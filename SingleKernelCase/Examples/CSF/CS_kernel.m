function f = CS_kernel(r)

% (c) XXXX

    
H=1;
beta=1/4;

f = H./(1 + r.^2).^beta;


end