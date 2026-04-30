% Authored by Jinchao Feng and Sui Tang
% smoke_test.m
% Tiny end-to-end check that the four paper experiments wire together:
% for each of {ODNS d=1, ODNS d=2, CSNS cut-off (paper CSNS1), CSNS smooth (paper CSNS2)}
% we (1) load the example, (2) generate ONE noise-free trajectory, (3) build a
% small (20-bin piecewise) feature matrix, (4) run a single FastLaplace fit
% with anchor index k=1, and (5) reconstruct phi_est at a few evaluation points.
% Total runtime ~30-90 s. No noise sweep, no model selection, no parfor.
%
% Pass criterion = pipeline completes without error and phi_est returns
% finite values at r > 0.

clear; close all;
addpaths;

cases = { ...
  struct('name','ODNS d=1 (step kernel)',         'kind','ODNS','d',1,'phi',@(r)ODNS_influence(r,1)); ...
  struct('name','ODNS d=2 (step kernel)',         'kind','ODNS','d',2,'phi',@(r)ODNS_influence(r,1)); ...
  struct('name','CSNS cut-off (paper CSNS1)',     'kind','CSNS','d',2,'phi',@(r)CS_kernel(r,2)); ...
  struct('name','CSNS smooth-decay (paper CSNS2)','kind','CSNS','d',2,'phi',@(r)CS_kernel(r,1)) };

passed = false(1,numel(cases));
times  = zeros(1,numel(cases));

for c = 1:numel(cases)
    C = cases{c};
    fprintf('\n========== smoke: %s ==========\n', C.name);
    t0 = tic;
    try
        %% (1) build sysInfo / obsInfo / solverInfo
        N = 50;                                    % half of paper's N to stay quick
        sysInfo.name      = C.kind;
        sysInfo.d         = C.d;
        sysInfo.N         = N;
        sysInfo.phi       = {C.phi};
        sysInfo.Noption   = 1;
        sysInfo.K         = 1;
        sysInfo.type_info = ones(N,1);
        sysInfo.type      = 1;
        sysInfo.RE        = [];
        sysInfo.flagxi    = 0;
        switch C.kind
            case 'ODNS'
                sysInfo.ode_order = 1;
                sysInfo.T_f       = 5;             % short horizon
                sysInfo.domain    = [0,10];
                sysInfo.mu0       = @() ODNS_init_config(C.d, N, C.d);
            case 'CSNS'
                sysInfo.ode_order = 2;
                sysInfo.T_f       = 5;
                sysInfo.domain    = [0,5];
                sysInfo.phi_type  = 'A';
                sysInfo.mu0       = @() CS_init_config(1,1, C.d, N, 1);
        end
        solverInfo.time_span = [0, sysInfo.T_f];
        solverInfo.option    = odeset('RelTol',1e-5,'AbsTol',1e-6);
        obsInfo.M               = 1;
        obsInfo.time_vec        = 0:1:sysInfo.T_f;
        obsInfo.use_derivative  = true;
        obsInfo.obs_noise       = 0;
        obsInfo.mu_trajnoise    = @(traj,sigma) trajUnifNoiseAdditive(traj,sigma);
        obsInfo.rho_T_histedges = linspace(0, sysInfo.domain(2), 200);
        fprintf('  setup OK\n');

        %% (2) generate one trajectory
        [~,~,dxpath_train, xpath_train] = Generate_training_data(sysInfo,obsInfo,solverInfo);
        L = size(xpath_train,2);
        fprintf('  generated 1 trajectory: x is %dx%dx1, dx is %dx%dx1\n', ...
                size(xpath_train,1), L, size(dxpath_train,1), L);

        %% (3) small piecewise basis
        nbasis = 20;
        rmax   = sysInfo.domain(2) - sysInfo.domain(1);
        edges  = linspace(0, rmax, nbasis+1);
        basis_funs = localBasisFn(edges, 0);
        fprintf('  built %d-bin piecewise basis on [0,%g]\n', nbasis, rmax);

        %% (4) assemble F (same convention as the base scripts) for the M=1, L snapshots
        d = sysInfo.d; order = sysInfo.ode_order;
        F_all = [];
        for l = 1:L
            for i = 1:N
                r_temp = xpath_train(1:d*N,l,1) - repmat(xpath_train((i-1)*d+1:i*d,l,1), N, 1);
                r_temp = reshape(r_temp, d, []);
                r_norm = sqrt(sum(r_temp.^2, 1));
                F_col  = zeros(nbasis, d);
                switch order
                    case 1   % ODNS: rhs uses dx - r_temp
                        rhs = reshape(repmat(dxpath_train((i-1)*d+1:i*d,l,1), N, 1) - r_temp(:), d, []);
                    case 2   % CSNS: rhs uses dvdot - v_temp
                        v_temp = xpath_train(d*N+1:d*N+d*N,l,1) - repmat(xpath_train(d*N+(i-1)*d+1:d*N+i*d,l,1), N, 1);
                        v_temp = reshape(v_temp, d, []);
                        rhs    = reshape(repmat(dxpath_train(d*N+(i-1)*d+1:d*N+i*d,l,1), N, 1) - v_temp(:), d, []);
                end
                for k = 1:nbasis
                    F_col(k,:) = sum(basis_funs{k}(repmat(r_norm, d, 1)) .* rhs, 2)';
                end
                F_all = [F_all F_col];
            end
        end
        F = F_all';
        F(isnan(F) | isinf(F)) = 0;
        fprintf('  assembled F: %dx%d (rows = d*N*L*M = %d, cols = K = %d)\n', size(F), d*N*L*1, nbasis);

        %% (5) one FastLaplace fit (anchor k* = 1)
        kstar = 1;
        eta   = F(:, kstar);
        Phi   = -F(:, [1:kstar-1, kstar+1:nbasis]);
        sigma2 = var(eta)/1e3 + eps;
        delta_La = 1e-10;
        [weights, used_ids] = FastLaplace(Phi, eta, sigma2, delta_La, [], []);
        fprintf('  FastLaplace returned %d/%d nonzero weights\n', numel(used_ids), nbasis-1);

        %% (6) reconstruct phi_est at a few r values
        ind = [1:kstar-1, kstar+1:nbasis];
        coef = zeros(nbasis,1); coef(kstar) = 1; coef(ind(used_ids)) = weights;
        phi_est = @(r) phi_piecewise_const(r, coef, rmax);
        rgrid = [0.1, 0.4, 0.6, 1.0, 2.0];
        vals  = phi_est(rgrid);
        if any(~isfinite(vals)), error('phi_est returned non-finite values'); end
        fprintf('  phi_est evaluated: '); fprintf('%.3g ', vals); fprintf('\n');

        passed(c) = true;
        times(c)  = toc(t0);
        fprintf('  PASS (%.2fs)\n', times(c));
    catch ME
        times(c) = toc(t0);
        fprintf('  FAIL after %.2fs: %s\n', times(c), ME.message);
    end
end

%% summary
fprintf('\n========== smoke-test summary ==========\n');
for c = 1:numel(cases)
    if passed(c), tag = 'PASS';
    else,         tag = 'FAIL'; end
    fprintf('  %-40s  %s  (%.2fs)\n', cases{c}.name, tag, times(c));
end
fprintf('\n%d / %d cases passed.  Total %.1fs.\n', sum(passed), numel(passed), sum(times));
if all(passed)
    fprintf('\nAll core pipelines (Examples + Generate_data + infer_kernel + FSBL_tools) are wired correctly.\n');
end
