function basis = construct_nonparametric_basis(R, num_basis_fun, basis_info)
% random feature basis 
% INPUT: number of basis function and information on type of random basis
% OUTPUT: a structure containing num_basis_fun many basis functions of given type 

n = num_basis_fun;
if length(R) < n
    error('cannot recoginize this type of random basis');
else
    randind = randsample(length(R),n);
end

if isfield(basis_info, 'sigma')
    sigma = basis_info.sigma;
else
    sigma = 1;
end

basis.r = R(randind);
basis.f = cell(1, n);
switch basis_info.kernel_type
    case 'gaussian'
        for ind = 1 : n
            basis.f{ind} = @(r) exp(-(R(randind(ind))-r).^2/(2*sigma^2));
        end
    
    case 'exp'
        for ind = 1 : n
            basis.f{ind} = @(r) exp(-abs(R(randind(ind))-r)/(2*sigma^2));
        end
    otherwise
        error('cannot recoginize this type of random basis');
end

end