function out = MaternKernel(x,y,omega,sigma,n,nu,partialOmega)

if nu == 1/2

    if ~partialOmega
        out = (1/n^2) * sigma^2 * exp(-1 * omega * abs(x - y));
    else
        out = -1 * (1/n^2) * sigma^2 * exp(-1 * omega * abs(x - y)) * (omega * abs(x - y));
    end

elseif nu == 3/2

    if ~partialOmega
        out = (1/n^2) * sigma^2 * (1 + sqrt(3) * abs(x - y) * omega) * exp(-1 * sqrt(3) * abs(x - y) * omega);
    else
        out = -1 * (1/n^2) * sigma^2 * exp(-1 * sqrt(3) * abs(x - y) * omega) * (3 * abs(x - y)^2 * omega^2);
    end

elseif nu == 5/2

    if ~partialOmega
        out = (1/n^2) * sigma^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (1 + sqrt(5)*abs(x - y)*omega + 5/3*abs(x - y)^2*omega^2);
    else
        out = -1 * (1/n^2) * sigma^2 * exp(-1 * omega * abs(x - y) * sqrt(5)) * (abs(x-y)^2*5*omega^2) * (1 + sqrt(5) * abs(x-y) * omega) / 3;
    end
    
elseif nu == 7/2

    if ~partialOmega
        out = (1/n^2) * sigma^2 * exp(-1 * omega * abs(x - y) * sqrt(7)) * (1 + sqrt(7)*abs(x - y)*omega + 2/5*7*abs(x - y)^2*omega^2 + 7/15 * sqrt(7) * abs(x-y)^3 * omega^3);
    else
        out = -1 * (1/n^2) * sigma^2 * exp(-1 * omega * abs(x - y) * sqrt(7)) * (7/5 * abs(x - y)^2 * omega^2) * (1 + sqrt(7) * abs(x - y) * omega + 7/3 * abs(x - y)^2 * omega^2);
    end

else

    fprintf('\n Error: This Matern parameter value is not yet fully implemented.');

end