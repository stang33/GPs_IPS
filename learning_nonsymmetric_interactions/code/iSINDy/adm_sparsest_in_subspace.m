function c = adm_sparsest_in_subspace(N, lambda, maxit, tol)
% Alternating Directions Method (ADM) for finding a sparse vector
% in the subspace span(N), following the implicit-SINDy philosophy.
%
% We solve approximately:
%   min ||q||_1
%   s.t. q = N z,  ||z||_2 = 1
%
% using alternating soft-thresholding / projection steps.
%
% INPUTS:
%   N      : orthonormal basis of approximate null space, size K x d
%   lambda : shrinkage threshold
%   maxit  : max iterations
%   tol    : stopping tolerance
%
% OUTPUT:
%   c      : sparse vector in span(N)

    [K, d] = size(N);

    if d == 0
        c = [];
        return;
    end

    % random initialization in the subspace coordinates
    z = randn(d,1);
    z = z / max(norm(z), eps);

    % initialize q and dual variable y
    q = N * z;
    y = zeros(K,1);

    beta = 1.0;

    for it = 1:maxit
        q_old = q;

        % q-update: soft threshold
        q = soft_threshold(N*z + y/beta, lambda/beta);

        % z-update: projection to unit sphere in subspace coordinates
        z_tilde = N' * (q - y/beta);
        nz = norm(z_tilde);
        if nz < 1e-14
            c = [];
            return;
        end
        z = z_tilde / nz;

        % dual update
        y = y + beta * (N*z - q);

        if norm(q - q_old) / max(norm(q_old), 1e-12) < tol
            break;
        end
    end

    c = N * z;
end

function y = soft_threshold(x, t)
    y = sign(x) .* max(abs(x) - t, 0);
end