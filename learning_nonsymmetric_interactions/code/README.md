# Reproducibility code for the JSC manuscript

> *A Sparse Bayesian Learning Algorithm for Estimation of Interaction Kernels in the Motsch--Tadmor Model.* — Jinchao Feng and Sui Tang.

This directory contains MATLAB code that reproduces every numerical example
in `main_SBL.tex` from synthetic data generated on the fly, including the
implicit-SINDy comparison. The package is self-contained.

## Layout

```
code/
├── README.md                                  # this file
├── addpaths.m                                 # adds helper folders to MATLAB path
│
├── Examples/                                  # system definitions and initial conditions
│   ├── CSNS_def.m, ODNS_def.m                 # paper sysInfo / obsInfo / solverInfo
│   ├── CSNS/{CS_kernel.m, CS_init_config.m}   # Cucker–Smale kernel and IC sampler
│   ├── ODNS/{ODNS_influence.m, ODNS_init_config.m}
│   ├── LoadExampleDefinitions.m, SelectExample.m
│   └── uniform_dist.m
│
├── Generate_data/                             # ODE integration + noise model
│   ├── Generate_training_data.m, Generate_rhoT.m
│   ├── RHSfn_NS.m, RHSfn_NS1.m,               # 1st-order normalized RHS (ODNS)
│   ├── RHSfn_2nd_NS.m, RHSfn_2nd_NS1.m, ...   # 2nd-order normalized RHS (CSNS)
│   ├── trajUnifNoise{Additive,Multiplicative}.m
│   ├── Exact_derivative.m, Pairwise_distance.m, rho_empirical.m, ...
│
├── basis_func/                                # all basis-function utilities
│   ├── localBasisFn.m                         # piecewise-constant basis
│   ├── basis_B_spline*.m                      # B-spline basis + derivatives
│   ├── construct_random_feature_basis.m       # RF (Gaussian) basis
│   ├── BasisMatrix_construction.m, Column_construction.m, RHSfn_basis.m, Recfns.m
│   ├── phi_piecewise_const.m                  # fast piecewise-constant kernel evaluator
│   └── mcos.m, mexp.m
│
├── FSBL_tools/                                # sparse-Bayesian solver (Laplace prior)
│   ├── FastLaplace.m                          # main: Babacan 2009 hierarchical Laplace
│   ├── ReweightedLaplace.m, test_W.m
│
├── SR_Tool/                                   # Tipping's SparseBayes (Gaussian flat prior)
│   └── SparseBayes.m + SB2_*.m
│
├── iSINDy/                                    # implicit-SINDy comparison (Mangan et al. 2016)
│   ├── run_implicit_sindy_homogeneous.m       # SVD null-space + ADM sparsification
│   ├── adm_sparsest_in_subspace.m
│   └── stlsq_solver.m
│
├── Plotting/                                  # trajectory + kernel visualization helpers
│
├── Run_Examples_Bayesian.m                    # **MAIN entry point** for kernel recovery + iSINDy comparison
│                                              #   (interactive: prompts for ODNS or CSNS)
│
├── Run_Examples_Bayesian_ODNS1_Msweep.m       # M=1..10 trajectory-error sweep, ODNS d=1
├── Run_Examples_Bayesian_ODNS2_Msweep.m       # M=1..10 trajectory-error sweep, ODNS d=2
├── Run_Examples_Bayesian_CSNS1_Msweep.m       # M=1..10 trajectory-error sweep, CSNS smooth (paper CSNS2)
├── Run_Examples_Bayesian_CSNS2_Msweep.m       # M=1..10 trajectory-error sweep, CSNS cut-off (paper CSNS1)
│
├── make_figure.m                              # ../images/{ODNSd1,ODNSd2,CSNS1,CSNS2}_traj_error.pdf
├── make_table.m                               # ../table_OD.tex, ../table_CS.tex
└── smoke_test.m                               # ~10s health check, runs all 4 examples once
```

> **Naming convention note.** In the paper, the second-order experiments are
> labeled *CSNS1 = cut-off kernel* (Example 1) and *CSNS2 = smooth-decay
> kernel* (Example 2). The on-disk `.mat` filenames in `../trajectory_data/`
> follow the historical convention where `CSNS1*.mat` was produced with
> `CS_kernel(r,1)` (smooth) and `CSNS2*.mat` with `CS_kernel(r,2)` (cut-off).
> The post-processing scripts (`make_figure.m`, `make_table.m`) handle this
> swap internally so that the *displayed* labels match the paper.

## Quick start

```matlab
% one-time path setup
cd <path-to>/JSCpaper/learning_nonsymmetric_interactions/code
addpaths
```

### Smoke test (~10 s)

```matlab
smoke_test
```
Expected: `4 / 4 cases passed.  Total ~10s.`. Validates Examples,
Generate_data, basis_func, FSBL_tools, iSINDy paths.

### Reproduce kernel-recovery figures (Figs 2, 3, 7--9, 13) and the iSINDy comparison

```matlab
Run_Examples_Bayesian
% prompts you to pick "ODNS" or "CSNS"; the chosen example's
% sysInfo (kernel, dimension) is loaded from Examples/{ODNS,CSNS}_def.m;
% then the script:
%   1. generates M training trajectories and adds noise,
%   2. assembles the implicit feature matrix F,
%   3. runs sparse-Bayesian model selection (FastLaplace + Gaussian SB),
%   4. runs implicit-SINDy comparison (run_iSINDy_compare = 1 by default),
%   5. saves and plots phi_est, phi_est_iS, plus uncertainty bands.
```

To switch CSNS variant (cut-off ↔ smooth) edit
`Examples/CSNS_def.m` line 11 (`@(r)CS_kernel(r,1)` smooth ↔
`@(r)CS_kernel(r,2)` cut-off). To switch ODNS dimension edit
`Examples/ODNS_def.m` line 9 (`sysInfo.d`).

### Reproduce trajectory-prediction-error figures (Figs 5, 12) and tables (Tabs 3-4, 6-7)

```matlab
Run_Examples_Bayesian_ODNS1_Msweep     % ~30 min on 10-core parpool
Run_Examples_Bayesian_ODNS2_Msweep     % ~50 min
Run_Examples_Bayesian_CSNS1_Msweep     % ~2 hours
Run_Examples_Bayesian_CSNS2_Msweep     % ~6 hours

make_table                              % writes ../table_OD.tex, ../table_CS.tex
make_figure                             % writes ../images/*_traj_error.pdf
```

Pre-computed `.mat` files for all `_Msweep` runs are already shipped in
`../trajectory_data/` (220 files, ~80 KB each). With those present,
`make_table` and `make_figure` complete in well under a minute and skip the
expensive sweeps.

## Paper figure / table → script map

| paper element                            | script(s) to run                                                   | output |
|------|------|------|
| Fig. ODtraj (1) trajectory profiles      | `Run_Examples_Bayesian` (no fit, just show traj)                   | inline |
| Fig. fig:OD_ex1 (2) ODNS kernel recovery | `Run_Examples_Bayesian` → ODNS                                     | `images/phi_ODNSN100d{1,2}M3L6sigma*_visualization.png` |
| Fig. fig:OD_errors (3)                   | sweep over M inside `Run_Examples_Bayesian` → ODNS                 | `images/ODNSd{1,2}_error.png` |
| Fig. fig:OD_runtime (4)                  | sweep over M                                                       | `images/ODNS{1,2}_time.png` |
| **Fig. fig:OD_traj_error (5)**           | `Run_Examples_Bayesian_ODNS{1,2}_Msweep.m` + `make_figure`         | `images/ODNSd{1,2}_traj_error.pdf` |
| Fig. CKtraj (6) CS trajectories          | `Run_Examples_Bayesian` → CSNS                                     | inline |
| Fig. fig:CK1Estimation (7) cut-off φ̂    | `Run_Examples_Bayesian` → CSNS with `CS_kernel(r,2)`               | `images/phi_CSNSN100d2M5L6sigma*_visualization.png` |
| Fig. fig:CK2Estimation (8) smooth φ̂     | `Run_Examples_Bayesian` → CSNS with `CS_kernel(r,1)`               | `images/phi_CSNS2N100d2M5L6sigma*_visualization.png` |
| Fig. fig:CS_errors (9)                   | sweep over M                                                       | `images/CSNS{1,2}_error.png` |
| Fig. fig:CK2_trajEstimation (10)         | trajectory plot inside CSNS run                                    | `images/CSNS2_trajs.png` |
| Fig. fig:CS_runtime (11)                 | sweep over M                                                       | `images/CSNS{1,2}_time.png` |
| **Fig. fig:CS_traj_error (12)**          | `Run_Examples_Bayesian_CSNS{1,2}_Msweep.m` + `make_figure`         | `images/CSNS{1,2}_traj_error.pdf` |
| Fig. fig:sindy_comparison (13)           | `Run_Examples_Bayesian` (with `run_iSINDy_compare = 1`, default)   | `images/phi_CSNS{1,2}*_comparison.png` |
| Tab. tab:OD1_rate (1) MS success rate    | sweep over noise inside `Run_Examples_Bayesian` → ODNS             | tabulated by hand |
| **Tab. tab:OD-d{1,2}-traj (3,4)**        | `Run_Examples_Bayesian_ODNS{1,2}_Msweep.m` + `make_table`          | `../table_OD.tex` |
| **Tab. tab:CS-{cutoff,smooth}-traj (6,7)** | `Run_Examples_Bayesian_CSNS{1,2}_Msweep.m` + `make_table`        | `../table_CS.tex` |

## Software environment

- MATLAB R2021b or newer (R2023a used in the published runs).
- Required toolboxes: Statistics & Machine Learning, Parallel Computing.
- Tested on Linux x86\_64.

## License

MIT-style permissive license, see `../LICENSE` (added on archival deposit).
