# ============================================================================
# SPHERE: Simulation Study
# ============================================================================
# This script demonstrates how to:
#   1. Generate synthetic spatial transcriptomics data under the SPHERE model
#   2. Fit the model using fit_sphere() and recover known parameters
#   3. Run replicated simulations for benchmarking
#
# The data-generating process mirrors the Stan model exactly:
#   Y_ij ~ Poisson(N_i * lambda_ij)
#   log(lambda_ij) = mu0 + Beta_j + epsilon_ij + eta_ij * I(Z_j = 2)
# where eta_ij is a GP spatial effect present only for SE genes.
# ============================================================================


# ----------------------------------------------------------------------------
# 0. Setup
# ----------------------------------------------------------------------------

library(SPHERE)
library(posterior)
library(bayesplot)
library(ggplot2)

theme_set(theme_bw())

#stan_model_path <- system.file("stan", "SPHERE_stan.stan", package = "SPHERE")
# Or use a local path:
 stan_model_path <- "SPHERE_stan.stan"


# ============================================================================
#
# PART A: SINGLE SIMULATION RUN
#
# ============================================================================

# ----------------------------------------------------------------------------
# A1. Define simulation parameters
# ----------------------------------------------------------------------------

num_spots <- 150
num_genes <- 50
prop      <- c(0.8, 0.2)     # 80% non-SE, 20% SE genes
G         <- 10               # number of pathway groups
seed      <- 124

# GP parameter vectors for SE genes (one value per SE gene)
n_se  <- round(num_genes * prop[2])
tau_g <- seq(0.5, by = 0.05, length.out = n_se)
ell_g <- seq(1.0, by = 0.03, length.out = n_se)

# Generate pathway group assignments and spatial coordinates
gene_grp <- generate_gene_grp(num_genes, G)
spots    <- generate_spatial_spots(
  num_spots, x1 = 2, x2 = 4, y1 = 5, y2 = 11, seed = 35
)


# ----------------------------------------------------------------------------
# A2. Generate synthetic data
# ----------------------------------------------------------------------------

sim_data <- generate_genedata_model(
  spots       = spots,
  num_genes   = num_genes,
  prop        = prop,
  G           = G,
  gene_grp    = gene_grp,
  rho         = 0.9,
  tau_beta    = 10,
  tau_gs      = tau_g,
  ell_gs      = ell_g,
  depth_model = "poisson",
  depth_lambda = 5,
  eps_sd      = 0.20,
  target_mean = 10,
  max_mean_target = 200,
  seed        = seed,
  verbose     = TRUE
)


# ----------------------------------------------------------------------------
# A3. Fit the SPHERE model
# ----------------------------------------------------------------------------

sim_fit <- fit_sphere(
  data_mat        = sim_data$Y,
  spot            = spots,
  gene_group      = gene_grp,
  stan_model_path = stan_model_path,
  iter_sampling   = 3500,
  iter_warmup     = 1500,
  chains          = 2,
  seed            = 1,
  knots           = 30
)


# ----------------------------------------------------------------------------
# A4. Evaluate parameter recovery
# ----------------------------------------------------------------------------

# Check convergence
cat("\nParameters with Rhat > 1.05:\n")
print(sim_fit$summary[sim_fit$summary$rhat > 1.05, ])

# Compare estimated Z to ground truth
z_rows   <- grep("^Z\\[", sim_fit$summary$variable)
z_est    <- round(sim_fit$summary$mean[z_rows])
z_true   <- sim_data$Z
accuracy <- mean(z_est == z_true)
cat("\nGene classification accuracy:", round(accuracy * 100, 1), "%\n")

# Save single-run results
save(sim_fit, sim_data, file = "SPHERE_simulation_single.RData")
cat("Results saved to SPHERE_simulation_single.RData\n")


# ============================================================================
#
# PART B: REPLICATED SIMULATION STUDY
#
# ============================================================================

# ----------------------------------------------------------------------------
# B1. Single-replication function
# ----------------------------------------------------------------------------

#' Run one simulation replication
#'
#' Generates data with a unique seed, fits SPHERE, and saves results.
#'
#' @param rep_id    Integer. Replication index.
#' @param num_spots Integer. Number of spatial spots.
#' @param num_genes Integer. Number of genes.
#' @param base_seed Integer. Base seed; actual seed = base_seed + rep_id.
#' @param out_dir   Character. Directory to save results.
#' @param stan_model_path Character. Path to the Stan model file.
#'
#' @return A list with fit results and ground-truth data (also saved to disk).
run_one_replication <- function(rep_id,
                                num_spots = 150,
                                num_genes = 50,
                                base_seed = 124,
                                out_dir   = "simulation_results",
                                stan_model_path = NULL) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  iter_seed <- base_seed + rep_id
  
  cat("\n========================================\n")
  cat("Replication", rep_id, "| Seed:", iter_seed, "\n")
  cat("========================================\n")
  
  # --- Simulation setup ---
  prop   <- c(0.8, 0.2)
  G      <- 10
  n_se   <- round(num_genes * prop[2])
  tau_g  <- seq(0.5, by = 0.05, length.out = n_se)
  ell_g  <- seq(1.0, by = 0.03, length.out = n_se)
  
  gene_grp <- generate_gene_grp(num_genes, G)
  spots    <- generate_spatial_spots(
    num_spots, x1 = 2, x2 = 4, y1 = 5, y2 = 11, seed = 35
  )
  
  # --- Generate data ---
  sim <- generate_genedata_model(
    spots = spots, num_genes = num_genes, prop = prop, G = G,
    gene_grp = gene_grp, rho = 0.9, tau_beta = 10,
    tau_gs = tau_g, ell_gs = ell_g,
    depth_model = "poisson", depth_lambda = 5,
    eps_sd = 0.20, target_mean = 10, max_mean_target = 200,
    seed = iter_seed, verbose = TRUE
  )
  
  # --- Fit model ---
  if (is.null(stan_model_path)) {
    stan_model_path <- system.file("stan", "SPHERE_stan.stan",
                                   package = "SPHERE")
  }
  
  fit_result <- fit_sphere(
    data_mat        = sim$Y,
    spot            = spots,
    gene_group      = gene_grp,
    stan_model_path = stan_model_path,
    iter_sampling   = 3500,
    iter_warmup     = 1500,
    chains          = 2,
    seed            = 1,
    knots           = 30
  )
  
  # --- Bundle and save ---
  output <- list(fit = fit_result, data = sim)
  
  fname <- sprintf("sim_rep%02d_n%d_p%d_seed%d.rds",
                   rep_id, num_spots, num_genes, iter_seed)
  saveRDS(output, file.path(out_dir, fname))
  cat("Saved:", fname, "\n")
  
  return(output)
}


# ----------------------------------------------------------------------------
# B2. Run replications
# ----------------------------------------------------------------------------

n_reps  <- 10
out_dir <- "simulation_results"

# Option 1: Sequential (good for debugging)
results <- lapply(seq_len(n_reps), function(i) {
  run_one_replication(
    rep_id    = i,
    num_spots = 150,
    num_genes = 50,
    base_seed = 124,
    out_dir   = out_dir
  )
})

# Option 2: Parallel (uncomment to use)
#
# library(future)
# library(future.apply)
# plan(multisession, workers = 4)
#
# results <- future_lapply(seq_len(n_reps), function(i) {
#   run_one_replication(
#     rep_id    = i,
#     num_spots = 150,
#     num_genes = 50,
#     base_seed = 124,
#     out_dir   = out_dir
#   )
# }, future.seed = TRUE)
#
# plan(sequential)


# ----------------------------------------------------------------------------
# B3. Aggregate results across replications
# ----------------------------------------------------------------------------

rds_files   <- list.files(out_dir, pattern = "^sim_rep.*\\.rds$",
                          full.names = TRUE)
all_results <- lapply(rds_files, readRDS)

cat("\nLoaded", length(all_results), "replications.\n")

# Classification accuracy per replication
accuracies <- sapply(all_results, function(res) {
  z_rows <- grep("^Z\\[", res$fit$summary$variable)
  z_est  <- round(res$fit$summary$mean[z_rows])
  z_true <- res$data$Z
  mean(z_est == z_true)
})

cat("\n--- Classification Accuracy Across Replications ---\n")
cat("Mean: ", round(mean(accuracies) * 100, 1), "%\n")
cat("SD:   ", round(sd(accuracies) * 100, 1), "%\n")
cat("Range:", round(min(accuracies) * 100, 1), "-",
    round(max(accuracies) * 100, 1), "%\n")