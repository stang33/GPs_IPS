function f = CS_kernel(r,type)

% Authored by Jinchao Feng and Sui Tang

switch type

    case 1
        H=1;
        beta=1/4;
        % beta=4;
        
        f = H./(1 + r.^2).^beta;

    case 2
        
        f = (r >= 0) & (r <= 0.5);
end    
end 