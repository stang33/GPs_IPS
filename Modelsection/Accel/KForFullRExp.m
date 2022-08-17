function out = KForFullRExp(x,y,omega,delta,n)

    out = (1/n^2) * delta^2 * exp(-1 * omega * abs(x - y));

end