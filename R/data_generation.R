## ============================================================
## SPHERE: Data Generation Functions
## ============================================================


## ------------------------------------------------------------
## Function to generate fixed number of spots (coordinates)
## ------------------------------------------------------------

#' Generate Spatial Spot Coordinates
#'
#' Generates a data frame of random spatial coordinates for a given number
#' of spots within a specified bounding box.
#'
#' @param num_spots Integer. Number of spatial spots to generate.
#' @param x1 Numeric. Minimum x-coordinate (default: 0).
#' @param x2 Numeric. Maximum x-coordinate (default: 1).
#' @param y1 Numeric. Minimum y-coordinate (default: 0).
#' @param y2 Numeric. Maximum y-coordinate (default: 1).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame with columns \code{x} and \code{y} containing
#'   the spatial coordinates of each spot.
#'
#' @importFrom stats runif
#' @export
generate_spatial_spots <- function(num_spots, x1 = 0, x2 = 1,
                                   y1 = 0, y2 = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data.frame(
    x = runif(num_spots, min = x1, max = x2),
    y = runif(num_spots, min = y1, max = y2)
  )
}



## ------------------------------------------------------------
## Function to simulate a named spatial pattern
## ------------------------------------------------------------

#' Simulate a Spatial Expression Pattern
#'
#' Generates a continuous spatial effect \eqn{\eta \in [0, 1]} for a set
#' of spatial spots using smooth spatial patterns. Unlike binary approaches,
#' this function returns a continuously varying spatial effect that more
#' realistically captures gradual expression changes across tissue.
#'
#' @param spots A data frame with columns \code{x} and \code{y} giving
#'   spatial coordinates.
#' @param pat_type Character. Type of spatial pattern. One of
#'   \code{"hotspot"}, \code{"streak"}, \code{"gradient"}, \code{"ring"},
#'   or \code{"wave"}. Use \code{NULL} for no pattern (returns zero eta).
#' @param high_exp_prop Numeric. Controls the spatial spread of the pattern —
#'   approximately the proportion of spots in the high-expression region
#'   (default: 0.3).
#' @param grad_percent Numeric. Gradient threshold in normalized coordinates
#'   for the \code{"gradient"} pattern (default: 0.3).
#' @param gene_id Integer or NULL. Gene index used to assign gene-specific
#'   random pattern centers. If NULL, uses the center of the tissue
#'   (default: NULL).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame with columns \code{x}, \code{y}, and \code{eta},
#'   where \code{eta} is a continuous spatial effect in \eqn{[0, 1]}.
#'
#' @importFrom stats runif
#' @export
simulate_expression <- function(spots, pat_type, high_exp_prop = 0.3,
                                grad_percent = 0.3, gene_id = NULL,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  x1 <- min(spots$x); x2 <- max(spots$x)
  y1 <- min(spots$y); y2 <- max(spots$y)
  x_range <- x2 - x1; y_range <- y2 - y1

  # Guard against zero range (single row/column of spots)
  if (x_range < 1e-12) x_range <- 1
  if (y_range < 1e-12) y_range <- 1

  # Normalize to [0, 1]
  x_norm <- (spots$x - x1) / x_range
  y_norm <- (spots$y - y1) / y_range
  n      <- nrow(spots)

  # Default: no spatial effect
  eta <- rep(0, n)
  if (is.null(pat_type)) {
    return(data.frame(x = spots$x, y = spots$y, eta = eta))
  }

  # Gene-specific random center
  if (!is.null(gene_id)) {
    cx <- runif(1, 0.15, 0.85)
    cy <- runif(1, 0.15, 0.85)
  } else {
    cx <- 0.5; cy <- 0.5
  }

  if (pat_type == "hotspot") {
    dist2 <- (x_norm - cx)^2 + (y_norm - cy)^2
    bw    <- sqrt(high_exp_prop / pi)
    eta   <- exp(-dist2 / (2 * bw^2))

  } else if (pat_type == "streak") {
    streak_x <- if (!is.null(gene_id)) runif(1, 0.15, 0.85) else cx
    bw  <- high_exp_prop / 2
    eta <- exp(-(x_norm - streak_x)^2 / (2 * bw^2))

  } else if (pat_type == "gradient") {
    base_mid  <- 1 - grad_percent
    midpoint  <- if (!is.null(gene_id)) runif(1, base_mid - 0.1,
                                              base_mid + 0.1) else base_mid
    steepness <- 10
    eta       <- 1 / (1 + exp(-steepness * (x_norm - midpoint)))

  } else if (pat_type == "ring") {
    dist2      <- (x_norm - cx)^2 + (y_norm - cy)^2
    target_r   <- sqrt(high_exp_prop / pi) * 1.2
    ring_width <- 0.08
    eta        <- exp(-(sqrt(dist2) - target_r)^2 / (2 * ring_width^2))

  } else if (pat_type == "wave") {
    phase       <- if (!is.null(gene_id)) runif(1, -1, 1) else 0
    freq        <- if (!is.null(gene_id)) runif(1, 2.5, 5.5) else 4
    wave_y      <- sin(freq * pi * (x_norm + phase))
    wave_scaled <- (wave_y - min(wave_y)) /
      (max(wave_y) - min(wave_y) + 1e-12)
    dist_to_curve <- abs(y_norm - wave_scaled)
    bw  <- 0.08
    eta <- exp(-dist_to_curve^2 / (2 * bw^2))

  } else {
    warning(paste("Unknown pat_type:", pat_type, "- returning zero eta."))
  }

  # Ensure eta is in [0, 1]
  eta <- pmax(0, pmin(1, eta))
  return(data.frame(x = spots$x, y = spots$y, eta = eta))
}


## ------------------------------------------------------------
## Function -- Pattern-based ST count data generator
## ------------------------------------------------------------

#' Generate Pattern-Based Spatial Transcriptomics Count Data
#'
#' Simulates spatial transcriptomics count data with a user-specified
#' spatial expression pattern. Genes are classified as spatially expressed
#' (SE) or non-spatially expressed (non-SE) according to \code{prop}.
#' SE genes follow a pattern defined by \code{pat_type}, while non-SE
#' genes have flat Poisson expression.
#'
#' @param spots A data frame or matrix of spatial coordinates
#'   (columns \code{x} and \code{y}).
#' @param num_genes Integer. Total number of genes to simulate.
#' @param prop Numeric vector of length 2 summing to 1. Proportions of
#'   non-SE and SE genes respectively (default: \code{c(0.8, 0.2)}).
#' @param pat_type Character. Spatial pattern type for SE genes. One of
#'   \code{"hotspot"}, \code{"streak"}, \code{"gradient"}, \code{"ring"},
#'   or \code{"wave"}.
#' @param seed Integer or NULL. Random seed for reproducibility
#'   (default: NULL).
#' @param G Integer. Number of biological pathway groups.
#' @param gene_grp Integer vector of length \code{num_genes} giving
#'   pathway membership for each gene.
#' @param cor_grp Numeric. CAR spatial correlation parameter \eqn{\rho}.
#' @param tau_beta Numeric. CAR precision parameter for pathway effects.
#' @param high_exp_prop Numeric. Proportion of spots in the
#'   high-expression region (default: 0.3).
#' @param grad_percent Numeric. Gradient threshold for the
#'   \code{"gradient"} pattern (default: 0.2).
#' @param beta0 Numeric. Log-scale baseline expression (default: 2).
#' @param beta1 Numeric or NULL. Log-scale high expression for SE genes.
#'   If NULL, defaults to \code{beta0 + log(3)} (3x fold change).
#' @param sigma Numeric. Standard deviation of log-scale noise
#'   (default: 0.3).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{Y}{Integer matrix (\eqn{n \times p}) of simulated counts.}
#'   \item{se_idx}{Integer vector of SE gene column indices.}
#'   \item{nonse_idx}{Integer vector of non-SE gene column indices.}
#'   \item{Beta}{Numeric vector of CAR pathway effects.}
#'   \item{epsilon}{Numeric matrix of log-scale noise.}
#'   \item{beta0}{Baseline log-mean used.}
#'   \item{beta1}{High-expression log-mean used.}
#' }
#'
#' @importFrom stats rnorm rpois quantile median
#' @export
generate_genedata_pattern <- function(
    spots, num_genes, prop = c(0.8, 0.2), pat_type, seed = NULL,
    G, gene_grp, cor_grp, tau_beta,
    high_exp_prop = 0.3,
    grad_percent  = 0.2,
    beta0         = 2,
    beta1         = NULL,
    sigma         = 0.3) {

  if (!is.null(seed)) set.seed(seed)

  if (length(prop) != 2 || any(prop < 0) || abs(sum(prop) - 1) > 1e-6)
    stop("prop must be a vector of length 2 summing to 1")

  spots <- as.data.frame(spots)
  colnames(spots)[1:2] <- c("x", "y")

  n <- nrow(spots)
  p <- num_genes

  if (is.null(beta1)) beta1 <- beta0 + log(3)

  num_cat    <- floor(p * prop)
  num_cat[which.max(prop * p - num_cat)] <-
    num_cat[which.max(prop * p - num_cat)] + (p - sum(num_cat))
  n_nonse <- num_cat[1]
  n_se    <- num_cat[2]

  se_idx    <- sort(sample(seq_len(p), n_se))
  nonse_idx <- setdiff(seq_len(p), se_idx)
  Z         <- rep(1L, p)
  Z[se_idx] <- 2L

  Beta    <- create_car_beta2(p, gene_grp, cor_grp, tau_beta)
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)

  result_mat <- matrix(0L, nrow = n, ncol = p)

  for (j in nonse_idx) {
    loglambda_j     <- beta0 + Beta[j] + epsilon[, j]
    result_mat[, j] <- rpois(n, lambda = exp(loglambda_j))
  }

  for (j in se_idx) {
    pattern_data <- simulate_expression(
      spots, pat_type = pat_type, high_exp_prop = high_exp_prop,
      grad_percent = grad_percent, gene_id = j)

    loglambda_j            <- beta0 + Beta[j] + epsilon[, j]
    high_mask              <- pattern_data$mark == "high"
    med_mask               <- pattern_data$mark == "medium"
    loglambda_j[high_mask] <- beta1 + Beta[j] + epsilon[high_mask, j]
    loglambda_j[med_mask]  <- (beta0 + beta1) / 2 + Beta[j] +
      epsilon[med_mask, j]

    result_mat[, j] <- rpois(n, lambda = exp(loglambda_j))
  }

  gene_names            <- character(p)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx]    <- paste0("Gene", se_idx, "_SE_", pat_type)
  colnames(result_mat)  <- gene_names

  cat("---- Pattern Simulation Summary ----\n")
  cat("n spots   :", n, " | P genes:", p, "\n")
  cat("SE genes  :", n_se, " | Non-SE genes:", n_nonse, "\n")
  cat("SE indices:", paste(se_idx, collapse = ", "), "\n")
  cat("Pattern   :", pat_type, "\n")
  cat("beta0 =", round(beta0, 3), "| beta1 =", round(beta1, 3),
      "| fold change =", round(exp(beta1 - beta0), 2), "x\n")
  cat("sigma     :", sigma, "\n")
  cat("Observed counts: median =", median(result_mat),
      "| 95% =", paste(round(quantile(as.vector(result_mat),
                                      c(.025, .975))), collapse = ", "),
      "| max =", max(result_mat), "\n")
  cat("------------------------------------\n")

  return(list(
    Y         = result_mat,
    se_idx    = se_idx,
    nonse_idx = nonse_idx,
    Beta      = Beta,
    epsilon   = epsilon,
    beta0     = beta0,
    beta1     = beta1
  ))
}


## ------------------------------------------------------------
## Function -- SPHERE model data generator
## ------------------------------------------------------------

#' Generate SPHERE Model Simulation Data
#'
#' Simulates spatial transcriptomics count data consistent with the
#' SPHERE model. Includes library size variation, CAR pathway effects,
#' gene-specific Gaussian Process spatial effects, and auto-calibrated
#' baseline expression.
#'
#' @param spots A matrix or data frame of spatial coordinates
#'   (at least 2 columns: x, y).
#' @param num_genes Integer. Total number of genes to simulate.
#' @param prop Numeric vector of length 2 summing to 1. Proportions of
#'   non-SE and SE genes respectively (default: \code{c(0.8, 0.2)}).
#' @param G Integer. Number of biological pathway groups.
#' @param gene_grp Integer vector of length \code{num_genes} giving
#'   pathway membership for each gene.
#' @param rho Numeric. CAR spatial correlation for pathway effects
#'   (default: 0.85).
#' @param tau_beta Numeric. CAR precision for pathway effects
#'   (default: 10).
#' @param tau_gs Numeric vector or NULL. GP amplitude for each SE gene.
#'   If NULL, sampled from \code{Uniform(0.05, 0.30)}.
#' @param ell_gs Numeric vector or NULL. GP lengthscale for each SE gene.
#'   If NULL, sampled from \code{Uniform(0.50, 2.00)}.
#' @param nugget Numeric. Nugget added to GP covariance diagonal for
#'   numerical stability (default: 1e-6).
#' @param eps_sd Numeric. Standard deviation of iid log-scale noise
#'   (default: 0.20).
#' @param depth_model Character. Model for sequencing depth. One of
#'   \code{"poisson"}, \code{"negbin"}, or \code{"fixed"}
#'   (default: \code{"poisson"}).
#' @param depth_lambda Numeric. Mean for Poisson or NegBin depth
#'   (default: 5).
#' @param depth_size Numeric. Overdispersion for NegBin depth —
#'   smaller values give more overdispersion (default: 5).
#' @param depth_fixed Numeric. Fixed library size when
#'   \code{depth_model = "fixed"} (default: 1).
#' @param target_mean Numeric. Desired median count for auto-calibrating
#'   the global intercept \code{mu0} (default: 10).
#' @param seed Integer or NULL. Random seed for reproducibility
#'   (default: NULL).
#' @param verbose Logical. If TRUE, prints a simulation summary
#'   (default: TRUE).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{Y}{Integer matrix (\eqn{n \times p}) of simulated counts.}
#'   \item{N_i}{Integer vector of per-spot library sizes.}
#'   \item{mu0}{Numeric. Auto-calibrated global log-expression intercept.}
#'   \item{Beta}{Numeric vector of CAR pathway effects.}
#'   \item{epsilon}{Numeric matrix of iid log-scale noise.}
#'   \item{eta}{Numeric matrix of GP spatial effects.}
#'   \item{Z}{Integer vector of true gene classifications (1=non-SE, 2=SE).}
#'   \item{se_idx}{Integer vector of SE gene indices.}
#'   \item{tau_gs}{Numeric vector of GP amplitudes for SE genes.}
#'   \item{ell_gs}{Numeric vector of GP lengthscales for SE genes.}
#'   \item{lambda}{Numeric matrix of latent Poisson rates.}
#'   \item{Y_mean}{Numeric matrix of expected counts.}
#' }
#'
#' @importFrom stats rnorm rpois rnbinom quantile median runif
#' @export
gen_genedata_model <- function(
    spots, num_genes, prop = c(0.8, 0.2), G, gene_grp,
    rho = 0.85, tau_beta = 10, tau_gs = NULL, ell_gs = NULL,
    nugget = 1e-6, eps_sd = 0.20,
    depth_model  = c("poisson", "negbin", "fixed"),
    depth_lambda = 5,
    depth_size   = 5,
    depth_fixed  = 1,
    target_mean  = 10,
    seed = NULL, verbose = TRUE) {

  depth_model <- match.arg(depth_model)

  stopifnot(is.matrix(spots) || is.data.frame(spots))
  spots <- as.matrix(spots)
  stopifnot(ncol(spots) >= 2)

  if (length(prop) != 2 || any(prop < 0) || abs(sum(prop) - 1) > 1e-8)
    stop("prop must be length-2 and sum to 1 (Non-SE, SE).")

  if (length(gene_grp) != num_genes)
    stop("gene_grp must have length num_genes.")

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(spots)
  P <- num_genes

  if (depth_model == "poisson") {
    N_i <- rpois(n, lambda = depth_lambda) + 1L
  } else if (depth_model == "negbin") {
    N_i <- rnbinom(n, mu = depth_lambda, size = depth_size) + 1L
  } else {
    N_i <- rep.int(as.integer(depth_fixed), n)
  }

  n_se      <- round(P * prop[2])
  se_idx    <- sort(sample(seq_len(P), n_se))
  nonse_idx <- setdiff(seq_len(P), se_idx)
  Z         <- rep(1L, P)
  Z[se_idx] <- 2L

  Beta    <- create_car_beta2(num_genes = P, gene_grp = gene_grp,
                              rho = rho, tau_beta = tau_beta)
  epsilon <- matrix(rnorm(n * P, mean = 0, sd = eps_sd), nrow = n, ncol = P)
  eta     <- matrix(0, nrow = n, ncol = P)

  if (n_se > 0) {
    if (is.null(tau_gs)) tau_gs <- runif(n_se, 0.05, 0.30)
    if (is.null(ell_gs)) ell_gs <- runif(n_se, 0.50, 2.00)

    if (length(tau_gs) != n_se) stop("tau_gs must have length = number of SE genes.")
    if (length(ell_gs) != n_se) stop("ell_gs must have length = number of SE genes.")

    eta_se        <- sim_gs_eta2(spots = spots[, 1:2, drop = FALSE],
                                 num_genes = n_se,
                                 tau_gs = tau_gs, ell_gs = ell_gs,
                                 nugget = nugget)
    eta[, se_idx] <- eta_se
  } else {
    tau_gs <- numeric(0)
    ell_gs <- numeric(0)
  }

  log_part <- sweep(epsilon + eta, 2, Beta, `+`)
  mu0      <- log(target_mean) - log(median(N_i)) - median(log_part)

  lambda <- exp(mu0 + log_part)
  Y_mean <- outer(N_i, rep(1, P)) * lambda
  Y      <- matrix(rpois(n * P, lambda = as.vector(Y_mean)),
                   nrow = n, ncol = P)

  gene_names            <- character(P)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx]    <- paste0("Gene", se_idx, "_SE_model")
  colnames(Y)           <- gene_names

  if (verbose) {
    cat("---- ST Simulation Summary ----\n")
    cat("n spots:", n, " | P genes:", P, "\n")
    cat("SE genes:", n_se, " (", round(100 * n_se / P, 1), "% )\n", sep = "")
    cat("Depth model:", depth_model,
        " | median N_i:", median(N_i),
        " | range:", paste(range(N_i), collapse = " - "), "\n")
    cat("mu0 (auto-calibrated):", signif(mu0, 4), "\n")
    cat("Expected mean counts : median =", signif(median(Y_mean), 4),
        " | 95% =", paste(signif(quantile(as.vector(Y_mean),
                                          c(.025, .975)), 4),
                          collapse = ", "),
        " | max =", signif(max(Y_mean), 4), "\n")
    cat("Observed counts      : median =", median(Y),
        " | 95% =", paste(quantile(as.vector(Y), c(.025, .975)),
                          collapse = ", "),
        " | max =", max(Y), "\n")
    cat("-------------------------------\n")
  }

  return(list(
    Y       = Y,
    spots   = spots |> as.data.frame(),
    N_i     = N_i,
    mu0     = mu0,
    Beta    = Beta,
    epsilon = epsilon,
    eta     = eta,
    Z       = Z,
    se_idx  = se_idx,
    tau_gs  = tau_gs,
    ell_gs  = ell_gs,
    lambda  = lambda,
    Y_mean  = Y_mean
  ))
}


## ------------------------------------------------------------
## Function to filter genes with less total counts
## ------------------------------------------------------------

#' Filter Genes with Low Total Counts
#'
#' Removes genes (columns) whose total count across all spots falls
#' below a specified threshold.
#'
#' @param df A matrix or data frame of gene expression counts.
#' @param threshold Numeric. Minimum total count required to keep a
#'   gene (default: 10).
#' @param keep_first Logical. If TRUE, always keeps the first column
#'   regardless of its sum (default: TRUE).
#'
#' @return A filtered \code{data.table} with low-count genes removed.
#'
#' @importFrom data.table as.data.table
#' @export
filter_columns_by_sum <- function(df, threshold = 10, keep_first = TRUE) {
  ..cols_keep <- NULL  # fix R CMD CHECK note
  dt <- as.data.table(df)
  if (keep_first) {
    numeric_part <- dt[, -1, with = FALSE]
    sums         <- colSums(numeric_part, na.rm = TRUE)
    cols_keep    <- c(names(dt)[1], names(numeric_part)[sums > threshold])
  } else {
    sums      <- colSums(dt, na.rm = TRUE)
    cols_keep <- names(dt)[sums > threshold]
  }
  return(dt[, ..cols_keep])
}
