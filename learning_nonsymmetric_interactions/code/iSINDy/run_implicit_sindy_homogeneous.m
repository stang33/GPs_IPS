function [c_best, info] = run_implicit_sindy_homogeneous(F_train, F_test, opts)
% Faithful implicit-SINDy style comparison:
%   find sparsest nonzero c such that F c ~= 0 is small
% by searching for sparse vectors in an approximate null space of F_train.
%
% INPUTS:
%   F_train : training implicit matrix
%   F_test  : testing implicit matrix
%   opts    : struct
%       .null_dims    : candidate null-space dimensions, e.g. 1:5
%       .lambda_grid  : shrinkage thresholds for ADM
%       .maxit
%       .tol
%       .verbose
%
% OUTPUTS:
%   c_best  : selected sparse coefficient vector
%   info    : diagnostic information

    if ~isfield(opts, 'null_dims'),   opts.null_dims = 1:5; end
    if ~isfield(opts, 'lambda_grid'), opts.lambda_grid = logspace(-4, 0, 20); end
    if ~isfield(opts, 'maxit'),       opts.maxit = 500; end
    if ~isfield(opts, 'tol'),         opts.tol = 1e-8; end
    if ~isfield(opts, 'verbose'),     opts.verbose = false; end

    [~, K] = size(F_train);

    % SVD of training matrix
    [~, S, V] = svd(F_train, 'econ');
    svals = diag(S);

    best_score = inf;
    c_best = zeros(K,1);

    info = struct();
    info.best_null_dim = NaN;
    info.best_lambda = NaN;
    info.best_residual = inf;
    info.best_score = inf;
    info.all = [];

    counter = 0;

    for dnull = opts.null_dims
        if dnull >= size(V,2)
            continue;
        end

        % approximate null space = span of right singular vectors
        % associated with the smallest singular values
        N = V(:, end-dnull+1:end);

        for lam = opts.lambda_grid
            counter = counter + 1;

            c = adm_sparsest_in_subspace(N, lam, opts.maxit, opts.tol);

            if isempty(c) || norm(c) < 1e-14
                continue;
            end

            % normalize for consistency
            c = c / max(norm(c), eps);

            % fix sign deterministically
            [~, idxm] = max(abs(c));
            if c(idxm) < 0
                c = -c;
            end

            % residual on test set
            resid_test = norm(F_test * c) / max(norm(c), eps);

            % sparsity
            nnz_c = nnz(abs(c) > 1e-8);

            % simple model-selection score: residual + small sparsity penalty
            score = resid_test + 1e-3 * nnz_c;

            info.all(counter).null_dim = dnull;
            info.all(counter).lambda = lam;
            info.all(counter).residual = resid_test;
            info.all(counter).nnz = nnz_c;
            info.all(counter).score = score;
            info.all(counter).c = c;

            if score < best_score
                best_score = score;
                c_best = c;
                info.best_null_dim = dnull;
                info.best_lambda = lam;
                info.best_residual = resid_test;
                info.best_score = score;
            end
        end
    end

    if opts.verbose
        fprintf('implicit-SINDy best score = %.4e\n', info.best_score);
    end
end