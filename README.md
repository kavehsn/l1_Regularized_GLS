# ℓ₁-Regularized GLS

[![DOI](https://img.shields.io/badge/DOI-10.1214%2F26--EJS2515-blue.svg)](https://doi.org/10.1214/26-EJS2515)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19211879.svg)](https://doi.org/10.5281/zenodo.19211879)

This repository accompanies the research paper:

> **ℓ₁-Regularized Generalized Least Squares**
> *Kaveh Salehzadeh-Nobari & Alex Gibberd (2026)*
> https://doi.org/10.1214/26-EJS2515

It provides all simulation code required to replicate the Monte Carlo experiments and figures in the paper.

---

## What This Repository Provides

- HPC-ready Monte Carlo simulation comparing LASSO, GLS-LASSO, and FGLS-LASSO estimators
- Estimation error analysis under AR(1) errors with ρ ∈ {0, 0.5, 0.9, 0.99}
- Sign and support recovery probability computation
- Verification of the R'R vs Σ_{MA(π)} spectral norm bound (Lemma 3)
- Covariance and Cholesky factor heatmaps for AR(1), AR(2), and AR(10) models
- Publication-ready figures matching the paper (Figures 1–4)

---

## Repository Structure

| File | Description |
|------|-------------|
| `HPC_estimation_error_simulation.R` | Main Monte Carlo simulation (Section 4, runs on HPC cluster) |
| `submit_estimation_error.sh` | PBS job submission script for HPC |
| `estimation_error_figures.R` | Generates Figure 3: estimation error plots (3×3 grid) |
| `sign_recovery_analysis.R` | Computes sign recovery probabilities (Figure 4) |
| `covariance_cholesky_heatmaps.R` | Generates Figure 1: Γ/v² and Ψ_q/v heatmaps |
| `RtR_MA_bound_verification.R` | Generates Figure 3 (appendix): R'R vs Σ_{MA} comparison |

---

## Getting Started

### 1️⃣ Clone the repository

```bash
git clone https://github.com/kavehsn/l1_Regularized_GLS.git
cd l1_Regularized_GLS
```

### 2️⃣ Download simulation data from Zenodo

The simulation output files (~2 GB) are hosted on Zenodo:

> **DOI:** [10.5281/zenodo.19211879](https://doi.org/10.5281/zenodo.19211879)

Download the `.RData` files and place them in a `Data` folder in your home directory:

```bash
mkdir -p ~/Data
mv ~/Downloads/l2ErrorTbl_*.RData ~/Data/
```

Your `~/Data/` folder should contain:

```
~/Data/
├── l2ErrorTbl_0.RData        # ρ = 0 (white noise baseline)
├── l2ErrorTbl_0.5.RData      # ρ = 0.5
├── l2ErrorTbl_0.9.RData      # ρ = 0.9
├── l2ErrorTbl_0.99.RData     # ρ = 0.99 (near unit root)
└── ...
```

### 3️⃣ Install R dependencies

```r
install.packages(c("glmnet", "MASS", "pracma", "ggplot2",
                    "latex2exp", "gridExtra", "fields",
                    "ltsa", "colorspace"))
```

### 4️⃣ Reproduce the figures

```r
# Figure 1: Covariance and Cholesky heatmaps
source("covariance_cholesky_heatmaps.R")

# Figure 3: Estimation error (3×3 grid)
source("estimation_error_figures.R")

# Figure 4: Sign recovery
source("sign_recovery_analysis.R")
```

---

## Regenerating Simulation Data (Optional)

If you wish to regenerate the simulation output from scratch rather than downloading from Zenodo, submit the HPC job:

```bash
qsub submit_estimation_error.sh
```

This runs `HPC_estimation_error_simulation.R` for ρ = 0. Edit the `Rho` argument in the script for other values of ρ. Each run produces ~250 MB of output and takes approximately 24–48 hours on a single core.

---

## Citation

If you use this code, please cite:

```bibtex
@article{salehzadeh2024l1gls,
  title={$\ell_1$-Regularized Generalized Least Squares},
  author={Salehzadeh-Nobari, Kaveh and Gibberd, Alex},
  journal={arXiv preprint arXiv:2405.10719},
  year={2024}
}
```

---

## License

This project is licensed under the MIT License.
