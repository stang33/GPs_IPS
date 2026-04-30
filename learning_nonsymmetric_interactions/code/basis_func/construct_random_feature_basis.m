function basis = construct_random_feature_basis(R, num_basis_fun, basis_info)
% random feature basis 
% INPUT: number of basis function and information on type of random basis
% OUTPUT: a structure containing num_basis_fun many basis functions of given type 

n = num_basis_fun;
if isfield(basis_info, 'mu')
    mu = basis_info.mu;
else
    mu = 0;
end
if isfield(basis_info, 'sigma')
    sigma = basis_info.sigma;
else
    sigma = 1;
end

basis.knots = linspace(0, R, 2);
basis.knotIdxs = ones(1,n,'uint32');
basis.f = cell(1, n);
omegas = normrnd(mu, sigma, [1,2*n]);
omegas_u = unique(omegas);
omegas = omegas_u(1:n);
basis.omegas = omegas;
% omegas = unifrnd(0, 20, [1,n]);

switch basis_info.feature_type
    case 'cos'
        ps = unifrnd(0, 2*pi, [1,n]);
        basis.ps = ps;
        for ind = 1 : n
            p = ps(ind);
            omega = omegas(ind);
            basis.f{ind} = @(r) mcos(r, omega, p);
        end
    case 'gaussian'
        if isfield(basis_info, 'normalize') && basis_info.normalize
            for ind = 1 : n
                omega = omegas(ind);
                omegasq = omega.^2;
                weight = 2.*abs(omega)/sqrt(pi);
%                 weight = 1;
%                 num = 7; % 7, 10
%                 weight = (2.^num).*(abs(omega).^(num+1))/(pi.^((num-1)/2));
                basis.f{ind} = @(r) mexp(r, omegasq, weight);
            end
        else
            for ind = 1 : n
                omega = omegas(ind);
                basis.f{ind} = @(r) mexp(r, omega, 1);
            end
        end
    otherwise
        error('cannot recoginize this type of random basis');
end
return