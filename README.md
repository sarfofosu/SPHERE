# SPHERE <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/sarfofosu/SPHERE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sarfofosu/SPHERE/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**SPHERE** (**S**patial **P**oisson **H**ierarchical mod**E**l with pathway-info**R**med g**E**ne networks) is an R package for identifying spatially expressed (SE) genes in spatial transcriptomics data.

SPHERE combines:
- A **Poisson log-normal likelihood** for modeling gene expression counts
- A **low-rank Gaussian Process** to capture spatial dependence across tissue locations
- A **pathway-informed Conditional Autoregressive (CAR) prior** to model biologically meaningful gene-gene dependencies
- A **two-component mixture model** for classifying genes as spatially expressed or non-spatially expressed

By jointly modeling spatial and gene-level dependencies informed by biological pathways, SPHERE improves both sensitivity and interpretability in detecting SE genes.

---

## Installation

SPHERE requires [R](https://www.r-project.org/) (>= 4.1.0) and [CmdStan](https://mc-stan.org/users/interfaces/cmdstan).

### Step 1 — Install CmdStan

```r
# Install cmdstanr
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# Install CmdStan (run once)
cmdstanr::install_cmdstan()
```

### Step 2 — Install SPHERE from GitHub

```r
# Install devtools if needed
install.packages("devtools")

# Install SPHERE
devtools::install_github("sarfofosu/SPHERE")
```

### Step 3 — Load the package

```r
library(SPHERE)
```

---

## Quick Start

```r
library(SPHERE)

# Fit the SPHERE model
result <- fit_sphere(
  data_mat   = counts_matrix,   # n x p matrix of gene expression counts
  spot       = coord_matrix,    # n x 2 matrix of spatial coordinates
  gene_group = pathway_labels,  # vector of pathway membership per gene
  chains     = 3,
  iter_sampling = 2000,
  iter_warmup   = 1000
)

# View posterior summary
head(result$summary)

# Posterior probability of each gene being spatially expressed
z_summary <- result$summary[grep("^Z\\[", result$summary$variable), ]
print(z_summary)
```

---

## Model Description

SPHERE models spatial transcriptomics count data using a hierarchical Bayesian framework:

```
Y_ij ~ Poisson(N_i * exp(lambda_ij))

lambda_ij = mu_0 + beta_j + f_j(s_i)

f_j(s_i) ~ GP(0, k_RBF)   [low-rank approximation with r knots]

beta_j ~ CAR(rho, tau)     [pathway-informed CAR prior]

Z_j ~ Mixture(pi_j)        [SE / non-SE classification]
```

Where:
- `Y_ij` is the count for gene `j` at spot `i`
- `N_i` is the library size at spot `i`
- `f_j` captures the spatial expression pattern for gene `j`
- `beta_j` is the gene-level effect informed by pathway membership
- `Z_j = 2` indicates gene `j` is spatially expressed

---

## Main Functions

| Function | Description |
|---|---|
| `fit_sphere()` | Fit the SPHERE model to spatial transcriptomics data |
| `plot_gene_data()` | Visualize spatial expression patterns |
| `filter_columns_by_sum()` | Filter low-count genes before model fitting |
| `make_rbf_basis()` | Construct low-rank RBF Gaussian Process basis |
| `make_rbf_dist()` | Compute knot distances for GP approximation |

---

## Citation

A manuscript describing SPHERE is currently in preparation. Please check back soon for a citation link.

> Sarfo Fosu, E. et al. (in preparation). SPHERE: Spatial Poisson Hierarchical modEl 
> with pathway-infoRmed gEne networks for spatial transcriptomics.

---

## License

This package is licensed under the [MIT License](LICENSE).

---

## Contact

For questions or issues, please open a [GitHub Issue](https://github.com/sarfofosu/SPHERE/issues) 
or contact Emmanuel Sarfo Fosu at emmanuelsarfo45@yahoo.com.
