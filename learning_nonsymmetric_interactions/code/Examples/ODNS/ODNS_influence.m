function influence = ODNS_influence(r, type)
% influence = OD_influence(r, type)

% Authored by Jinchao Feng and Sui Tang

influence = zeros(size(r));
cutoff    = 0.5;
support   = 1;

switch type
    case 1
        ind            = 0 <= r & r <cutoff;
        influence(ind) = 1;
        ind            = cutoff <= r & r < support;
        influence(ind) = 0.1;
    case 2
        influence       = exp(-2*r)+exp(-4*r);
        %influence       = sin(r).^2;
        
    case 3
        influence       = 40*((-12)*r.^(-14) - (-6)*r.^(-8));
        %influence       = -4*r.^(-5) + 4*r.^(-3);

    case 4
        influence       = (r./10).^3 + (r./10).^5;

    case 5
        influence       = cos(r).^2 + cos(3*r).^2;   

    case 6
        influence       = - r.^(1) + r.^3;   

    case 10
        influence       = (1+r.^2).^(-1/4);     
    case 11
        influence       = (1+r.^2).^(-4);

    otherwise
end

return