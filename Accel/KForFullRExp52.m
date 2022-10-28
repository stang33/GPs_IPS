function out = KForFullRExp52(x,y,omega,delta,n)

    out = (1/n^2) * delta^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (1 + sqrt(5)*abs(x - y)*omega + 5/3*abs(x - y)^2*omega^2);

end