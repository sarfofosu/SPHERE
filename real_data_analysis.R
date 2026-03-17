# ============================================================================
# SPHERE: Spatial Variable Gene Detection in Human Breast Cancer
# ============================================================================
# This script demonstrates a complete SPHERE analysis workflow using
# human breast cancer spatial transcriptomics data. It covers:
#   1. Loading and preparing the expression data
#   2. Handling missing pathway annotations
#   3. Filtering pathways by size
#   4. Fitting the SPHERE model
#   5. Inspecting results
# ============================================================================


# ----------------------------------------
# 0. Load required packages
# ----------------------------------------

if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(
  purrr, furrr, knitr, readr, parallel, rstan, dplyr, cmdstanr, posterior,
  bayesplot, ggplot2, magrittr, stringr, tibble, tidyr, SPARK, HEARTSVG)


# ----------------------------------------
# 1. Load data
# ----------------------------------------

# Full human breast cancer spatial transcriptomics dataset
load("humanBC_dataset.RData")
# Pre-processed subset with selected genes and spots
load("BrC_New_set.RData")


# ----------------------------------------
# 2. Prepare expression data and spatial coordinates
# ----------------------------------------

# Expression count matrix (spots x genes)
# round() ensures integer counts for the Poisson likelihood
data_mat <- round(BrC_final)
# Spatial coordinates for each spot (x, y)
spot <- BrC_spot
# Quick sanity check
cat("Spots:", nrow(data_mat), "\n")
cat("Genes:", ncol(data_mat), "\n")



# ----------------------------------------
# 3. Prepare pathway annotations
# ----------------------------------------

# pathway_df is loaded from the .RData file and must contain a column
# named "Pathway" with one entry per gene (same order as columns of data_mat)

# 3a. Handle genes with missing pathway annotations
#     In the Stan CAR prior, a gene with no neighbors (n_j = 0) receives
#     an independent Normal(0, 10) prior on Beta[j], so these genes are
#     effectively treated as isolated.
na_idx <- which(is.na(pathway_df$Pathway))

if (length(na_idx) > 0) {
  # Give each unannotated gene a unique pathway label: "isolated_1", "isolated_2", ...
  pathway_df$Pathway[na_idx] <- paste0("isolated_", seq_along(na_idx))
  cat("Assigned", length(na_idx),
      "unannotated genes to individual singleton pathways.\n")
}
# 3b. Filter pathways by size 
gen_path <- filter_pathways_by_limit(pathway_df)
# Final gene-level pathway vector (length = number of genes)
gene_pathway_final <- gen_path$Pathway
cat("Pathways retained:", length(unique(gene_pathway_final)), "\n")


# ----------------------------------------
# 4. Specify Stan model path
# ----------------------------------------

# If running from a local file instead
stan_model_path <- "SPHERE_stan.stan"

# Point to the compiled Stan file shipped with the package (or a local copy)
#stan_model_path <- system.file("stan", "SPHERE_stan.stan", package = "SPHERE")


# ----------------------------------------
# 5. Fit the SPHERE model
# ----------------------------------------

#NB Runing this may take hours to complete
BrC_fit <- fit_sphere(
  data_mat       = data_mat,
  spot           = spot,
  gene_group     = gene_pathway_final,
  stan_model_path = stan_model_path,
  iter_sampling  = 5000,          # posterior draws per chain
  iter_warmup    = 2000,          # warmup (burn-in) per chain
  chains         = 3,             # number of independent chains
  seed           = 8,             # for reproducibility
  knots          = 30             # inducing points for low-rank GP
)

cat("Model fitting completed in", round(BrC_fit$runtime, 1), "seconds.\n")

# ----------------------------------------
# 6. Inspect results
# ----------------------------------------

# 6a. Check convergence diagnostics (Rhat should be < 1.05)
summary_df <- BrC_fit$summary
cat("\nParameters with Rhat > 1.05:\n")
print(summary_df[summary_df$rhat > 1.05, ])

# 6b. Extract posterior classification of genes
#     Z[j] = 1 means gene j is non-se Z[j] = 2 means gene j is se
z_rows <- grep("^Z\\[", summary_df$variable)
z_summary <- summary_df[z_rows, ]

# 6c. Save results
save(BrC_fit, z_summary, file = "SPHERE_BrC_results.RData")
cat("\nResults saved to SPHERE_BrC_results.RData\n")