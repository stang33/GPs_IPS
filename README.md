# GPs_IPS — Gaussian-Process and Sparse-Bayesian methods for Interacting Particle Systems

This repository collects the code and reproducibility packages accompanying our work on data-driven discovery of interaction kernels in interacting particle systems (IPS).

## Contents

### 1. Symmetric IPS — Gaussian-process kernel learning *(repo root)*

Original Gaussian-process framework for learning *symmetric* interaction kernels:

- `Run_GP_Examples.m` — synthetic-system entry point
- `Run_GP_RealData.m` — fish-school real-data example
- Helper folders: `infer_kernel/`, `Modelsection/`, `Plotting/`, `SingleKernelCase/`

### 2. Non-symmetric (Motsch–Tadmor) IPS — Sparse Bayesian Learning

Folder: [`learning_nonsymmetric_interactions/`](learning_nonsymmetric_interactions/)

Companion code and full reproducibility package for

> **A Sparse Bayesian Learning Algorithm for Estimation of Interaction Kernels in the Motsch–Tadmor Model** — Jinchao Feng and Sui Tang, submitted to the *Journal of Scientific Computing*.

Contents:
- `main_SBL.tex` / `main_SBL.pdf` — manuscript (23 pp)
- `response_letter.tex` / `response_letter.pdf` — response to reviewers
- `cover letter.tex`
- `table_OD.tex`, `table_CS.tex` — auto-generated trajectory-error tables
- `reference.bib`, Springer JSC class files (`spmpsci.bst`, `svjour3.cls`, `svglov3.clo`)
- `images/` — 35 figures used by the manuscript
- `code/` — self-contained MATLAB reproducibility code (see [`learning_nonsymmetric_interactions/code/README.md`](learning_nonsymmetric_interactions/code/README.md))
- `trajectory_data/` — 220 `.mat` files of pre-computed M-sweep results

### Quick reproduction (Sec. 2)

```bash
cd learning_nonsymmetric_interactions/code
matlab -batch "addpaths; smoke_test"            # ~10 s health check
matlab -batch "Params.ExampleName='ODNS'; Run_Examples_Bayesian"   # paper §4.1
matlab -batch "Params.ExampleName='CSNS'; Run_Examples_Bayesian"   # paper §4.2 + iSINDy comparison
```

## Citation

```bibtex
@unpublished{feng_tang_motsch_tadmor_2026,
  author = {Feng, Jinchao and Tang, Sui},
  title  = {A Sparse Bayesian Learning Algorithm for Estimation of Interaction Kernels in the Motsch--Tadmor Model},
  note   = {Submitted to the Journal of Scientific Computing},
  year   = {2026}
}
```

## License

Author-written MATLAB files are released under an MIT-style permissive license. Third-party modules retain their original copyright statements (Tipping for `code/SR_Tool/`, Babacan for `code/FSBL_tools/FastLaplace.m`, Mangan et al. for the iSINDy basis used in `code/iSINDy/`).
