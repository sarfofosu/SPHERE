# ============================================================================
# SPHERE: Simulation Study Example
# ============================================================================
# This script demonstrates how to:
#   1. Generate synthetic spatial transcriptomics data under the SPHERE model
#   2. Fit the SPHERE model using fit_sphere()
#   3. Evaluate parameter recovery and gene classification
#   4. Run replicated simulations for benchmarking
#
# The data-generating process mirrors the Stan model exactly:
#   Y_ij ~ Poisson(N_i * lambda_ij)
#   log(lambda_ij) = mu0 + Beta_j + epsilon_ij + eta_ij * I(Z_j = 2)
# where eta_ij is a GP spatial effect present only for SE genes.
# ============================================================================


# ----------------------------------------------------------------------------
# 0. Load package
# ----------------------------------------------------------------------------
library(SPHERE)


# ============================================================================
# PART A: SINGLE SIMULATION RUN
# ============================================================================

# ----------------------------------------------------------------------------
# A1. Define simulation parameters
# ----------------------------------------------------------------------------

num_spots <- 50
num_genes <- 20
prop      <- c(0.8, 0.2)     # 80% non-SE, 20% SE genes
G         <- 10               # number of pathway groups
seed      <- 124

# GP parameter vectors for SE genes (one value per SE gene)
n_se  <- round(num_genes * prop[2])
t_gs  <- seq(1.5, by = 0.05, length.out = n_se)
l_gs  <- seq(1.0, by = 0.03, length.out = n_se)

# Generate pathway group assignments and spatial coordinates
gene_grp <- generate_gene_grp(num_genes, G)
spots    <- generate_spatial_spots(
  num_spots, x1 = 2, x2 = 4, y1 = 5, y2 = 11, seed = 35
)


# ----------------------------------------------------------------------------
# A2. Generate synthetic data
# ----------------------------------------------------------------------------

sim_data <- generate_genedata_model(
  spots        = spots,
  num_genes    = num_genes,
  prop         = prop,
  G            = G,
  gene_grp     = gene_grp,
  rho          = 0.9,
  tau_beta     = 20,
  tau_gs       = t_gs,
  ell_gs       = l_gs,
  depth_model  = "negbin",
  depth_lambda = 5,
  depth_size   = 5,
  eps_sd       = 0.20,
  target_mean  = 10,
  seed         = seed,
  verbose      = TRUE
)

# Extract count matrix
count_mat <- sim_data$Y
cat("Count matrix dimensions:", nrow(count_mat), "spots x",
    ncol(count_mat), "genes\n")
cat("True SE genes:", paste(colnames(count_mat)[sim_data$se_idx],
                            collapse = ", "), "\n")


# ----------------------------------------------------------------------------
# A3. Fit the SPHERE model
# ----------------------------------------------------------------------------

sim_fit <- fit_sphere(
  data_mat          = count_mat,
  spot              = spots,
  gene_group        = gene_grp,
  iter_sampling     = 200,
  iter_warmup       = 200,
  chains            = 3,
  seed              = 8,
  knots             = 30,
  alpha             = c(10, 3),
  refresh           = 50,
  # Observation noise prior
  mu_noise          = 0,
  sd_noise          = 1,
  # GP lengthscale prior
  mu_gp_lengthscale = 0,
  sd_gp_lengthscale = 3,
  # GP amplitude prior
  mu_gp_amplitude   = 0,
  sd_gp_amplitude   = 12,
  # Global intercept prior
  mu_intercept      = 0,
  sd_intercept      = 1,
  # CAR correlation prior
  shape_beta_rho    = 5,
  rate_beta_rho     = 2,
  # CAR precision prior
  mu_beta_sig       = 0,
  sd_beta_sig       = 1
)


# ----------------------------------------------------------------------------
# A4. Evaluate results
# ----------------------------------------------------------------------------

# --- Convergence diagnostics ---
cat("\nParameters with Rhat > 1.05:\n")
high_rhat <- sim_fit$summary[!is.na(sim_fit$summary$rhat) &
                               sim_fit$summary$rhat > 1.05, ]
if (nrow(high_rhat) == 0) {
  cat("None — all parameters converged!\n")
} else {
  print(high_rhat)
}

# --- Gene classification results ---
z_rows <- grep("^Z\\[", sim_fit$summary$variable)
z_summary <- sim_fit$summary[z_rows, c("variable", "mean", "median")]
cat("\nPosterior Z summary (prob of being SE):\n")
print(z_summary)

# --- Classification accuracy ---
z_est    <- round(sim_fit$summary$mean[z_rows])
z_true   <- sim_data$Z
accuracy <- mean(z_est == z_true)
cat("\nGene classification accuracy:", round(accuracy * 100, 1), "%\n")

# --- SE gene detection ---
cat("\nTrue SE gene indices   :", sim_data$se_idx, "\n")
cat("Estimated SE gene indices:", which(z_est == 2), "\n")

# --- Runtime ---
cat("\nTotal runtime:", round(sim_fit$runtime, 1), "seconds\n")


# ============================================================================
# PART B: REPLICATED SIMULATION STUDY
# ============================================================================

# ----------------------------------------------------------------------------
# B1. Single replication function
# ----------------------------------------------------------------------------

#' Run One Simulation Replication
#'
#' Generates synthetic ST data with a unique seed, fits SPHERE,
#' evaluates classification accuracy, and saves results to disk.
#'
#' @param rep_id    Integer. Replication index.
#' @param num_spots Integer. Number of spatial spots (default: 50).
#' @param num_genes Integer. Number of genes (default: 20).
#' @param base_seed Integer. Base seed; actual seed = base_seed + rep_id.
#' @param out_dir   Character. Directory to save results.
#'
#' @return A list with \code{fit} and \code{data} elements.
run_one_replication <- function(rep_id,
                                num_spots = 50,
                                num_genes = 20,
                                base_seed = 124,
                                out_dir   = "simulation_results") {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  iter_seed <- base_seed + rep_id

  cat("\n========================================\n")
  cat("Replication", rep_id, "| Seed:", iter_seed, "\n")
  cat("========================================\n")

  # Simulation setup
  prop  <- c(0.8, 0.2)
  G     <- 10
  n_se  <- round(num_genes * prop[2])
  t_gs  <- seq(1.5, by = 0.05, length.out = n_se)
  l_gs  <- seq(1.0, by = 0.03, length.out = n_se)

  gene_grp <- generate_gene_grp(num_genes, G)
  spots    <- generate_spatial_spots(
    num_spots, x1 = 2, x2 = 4, y1 = 5, y2 = 11, seed = 35
  )

  # Generate data
  sim <- gen_genedata_model(
    spots        = spots,
    num_genes    = num_genes,
    prop         = prop,
    G            = G,
    gene_grp     = gene_grp,
    rho          = 0.9,
    tau_beta     = 20,
    tau_gs       = t_gs,
    ell_gs       = l_gs,
    depth_model  = "negbin",
    depth_lambda = 5,
    depth_size   = 5,
    eps_sd       = 0.20,
    target_mean  = 10,
    seed         = iter_seed,
    verbose      = FALSE
  )

  # Fit SPHERE model
  fit_result <- fit_sphere(
    data_mat          = sim$Y,
    spot              = spots,
    gene_group        = gene_grp,
    iter_sampling     = 200,
    iter_warmup       = 200,
    chains            = 3,
    seed              = 1,
    knots             = 30,
    alpha             = c(10, 3),
    refresh           = 0,
    mu_noise          = 0,   sd_noise          = 1,
    mu_gp_lengthscale = 0,   sd_gp_lengthscale = 3,
    mu_gp_amplitude   = 0,   sd_gp_amplitude   = 12,
    mu_intercept      = 0,   sd_intercept      = 1,
    shape_beta_rho    = 5,   rate_beta_rho     = 2,
    mu_beta_sig       = 0,   sd_beta_sig       = 1
  )

  # Evaluate classification accuracy
  z_rows   <- grep("^Z\\[", fit_result$summary$variable)
  z_est    <- round(fit_result$summary$mean[z_rows])
  accuracy <- mean(z_est == sim$Z)
  cat("Classification accuracy:", round(accuracy * 100, 1), "%\n")
  cat("Runtime:", round(fit_result$runtime, 1), "seconds\n")

  # Bundle and save
  output <- list(fit = fit_result, data = sim, accuracy = accuracy)

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
    num_spots = 50,
    num_genes = 20,
    base_seed = 124,
    out_dir   = out_dir
  )
})


# Option 2: Parallel (uncomment to use)
# library(future)
# library(future.apply)
# plan(multisession, workers = 4)
#
# results <- future_lapply(seq_len(n_reps), function(i) {
#   run_one_replication(
#     rep_id    = i,
#     num_spots = 50,
#     num_genes = 20,
#     base_seed = 124,
#     out_dir   = out_dir
#   )
# }, future.seed = TRUE)


# ----------------------------------------------------------------------------
# B3. Summarize replication results
# ----------------------------------------------------------------------------

# Extract accuracy across replications
accuracies <- sapply(results, function(r) r$accuracy)
runtimes   <- sapply(results, function(r) r$fit$runtime)

cat("\n===== Simulation Study Summary =====\n")
cat("Replications         :", n_reps, "\n")
cat("Mean accuracy        :", round(mean(accuracies) * 100, 1), "%\n")
cat("SD accuracy          :", round(sd(accuracies) * 100, 1), "%\n")
cat("Min accuracy         :", round(min(accuracies) * 100, 1), "%\n")
cat("Max accuracy         :", round(max(accuracies) * 100, 1), "%\n")
cat("Mean runtime (secs)  :", round(mean(runtimes), 1), "\n")
cat("=====================================\n")
