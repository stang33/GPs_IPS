function f = FM_kernel(r,N)

% (c) XXXX
    
Ca = 4;
la = 4;
Cr = 0.5;%0.36;%0.5;%C=Cr/Ca=0.9
lr = 0.5;%0.8;%0.5; %l=lr/la=0.5;

r0=0.05;
f       = zeros(size(r));
ind     = r>=r0;
f(ind) = N./r(ind).*(Ca/la*exp(-r(ind)/la) - Cr/lr*exp(-r(ind)/lr));

f0= N./(r0)*(Ca/la*exp(-r0/la)-Cr/lr*exp(-r0/lr));
f1= N./r0.*(Ca/la*exp(-r0/la)*(-1/la) + Cr/lr*exp(-r0/lr)/lr);
f2= -f0/r0;

ind  =  r<r0;
a=-((f1+f2)/f0)/(12*r0^11);
b=f0*exp(a*r0.^12);
f(ind)=b.*exp(-a*r(ind).^12);
end