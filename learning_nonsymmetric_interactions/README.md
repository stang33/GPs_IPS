# Reproducibility code: SBL for asymmetric interaction kernels

Companion code release for

> **A Sparse Bayesian Learning Algorithm for Estimation of Interaction Kernels in the Motsch–Tadmor Model** — Jinchao Feng and Sui Tang, submitted to the *Journal of Scientific Computing* (2026).

This folder contains a self-contained MATLAB reproducibility package that regenerates every numerical figure and table reported in the manuscript from synthetic data. The manuscript, response letter, and other paper-level materials are not redistributed here while the paper is under review; once published, a link to the official version will be added.

## Layout

```
learning_nonsymmetric_interactions/
├── README.md                       this file
├── code/                           reproducibility scripts (~80 .m files, see code/README.md)
└── trajectory_data/                220 .mat files of pre-computed M-sweep results
                                    (lets `make_table` and `make_figure` finish in seconds
                                     without re-running the multi-hour SBR sweeps)
```

## Quick start

```matlab
cd learning_nonsymmetric_interactions/code
addpaths
smoke_test                                            % ~10 s health check
Params.ExampleName = 'ODNS'; Run_Examples_Bayesian   % paper §4.1
Params.ExampleName = 'CSNS'; Run_Examples_Bayesian   % paper §4.2 + iSINDy comparison
Run_Examples_Bayesian_ODNS1_Msweep                    % trajectory error sweep, ODNS d=1
Run_Examples_Bayesian_ODNS2_Msweep                    %                      , ODNS d=2
Run_Examples_Bayesian_CSNS1_Msweep                    %                      , CSNS smooth
Run_Examples_Bayesian_CSNS2_Msweep                    %                      , CSNS cut-off
make_table                                            % regenerate trajectory-error tables
make_figure                                           % regenerate trajectory-error figures
```

See [`code/README.md`](code/README.md) for a per-figure / per-table reproduction recipe.

## Software environment

* MATLAB R2021b or newer (R2023a used in the published runs).
* Required toolboxes: Statistics & Machine Learning, Parallel Computing.
* Tested on Linux x86_64.

## Citation

Once the paper is published the citation will be updated; for now please cite the preprint:

```bibtex
@unpublished{feng_tang_motsch_tadmor_2026,
  author = {Feng, Jinchao and Tang, Sui},
  title  = {A Sparse {B}ayesian Learning Algorithm for Estimation of
            Interaction Kernels in the {M}otsch--{T}admor Model},
  note   = {Submitted to the Journal of Scientific Computing},
  year   = {2026}
}
```

## Code authorship

Author-written MATLAB files carry the header `% Authored by Jinchao Feng and Sui Tang`. Third-party modules retain their original copyright statements:

* `code/SR_Tool/`        — Mike Tipping (`SparseBayes`, Vector Anomaly Ltd, 2009)
* `code/FSBL_tools/FastLaplace.m` — S. Derin Babacan et al. (Northwestern, 2008)
* `code/iSINDy/`         — based on Mangan et al., *J. Roy. Soc. Interface* (2016); ADM null-space implementation contributed by Jinchao Feng
