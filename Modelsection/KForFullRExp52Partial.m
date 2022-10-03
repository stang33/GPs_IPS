function out = KForFullRExp52Partial(x,y,omega,delta,n)

    %out = (1/n^2) * delta^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (1 + sqrt(5)*abs(x - y)*omega + 5/3*abs(x - y)^2*omega^2);

    %out = (-1 * abs(x - y) * sqrt(5)) * (1/n^2) * delta^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (1 + sqrt(5)*abs(x - y)*omega + 5/3*abs(x - y)^2*omega^2) + (sqrt(5)*abs(x - y) + 10/3 * abs(x - y)^2 * omega) * (delta^2 / n^2) * exp(-1 * omega * abs(x-y) * sqrt(5));

    %out = (delta/n^2) * exp(-1 * sqrt(5) * abs(x-y) * omega) *(-5/3 * abs(x-y)^2 * omega - (5* sqrt(5)) / 3 * abs(x-y)^3 * omega^2);
    %out = (1/n^2) * delta^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (1 + sqrt(5)*abs(x - y)*omega + 5/3*abs(x - y)^2*omega^2); %orig from kphimatern
    out = -1 * (1/n^2) * delta^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (abs(x-y)^2*5*omega^2) * (1 + sqrt(5) * abs(x-y) * omega) / 3;%(1 + sqrt(5)*abs(x - y)*omega + 5/3*abs(x - y)^2*omega^2);

end