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
#' Simulates spatial transcriptomics (ST) count data with configurable
#' spatial expression patterns, sequencing depth variation, pathway-level
#' dependence, and realistic count calibration. Genes are classified as
#' spatially expressed (SE) or non-spatially expressed (non-SE) according
#' to \code{prop}. SE genes follow spatial patterns defined by
#' \code{pat_type}, while non-SE genes exhibit baseline expression with
#' controlled zero inflation.
#'
#' This generator incorporates:
#' \itemize{
#'   \item Pathway-level correlation via Conditional Autoregressive (CAR) effects
#'   \item Spatial signal scaling through \code{eta_scale}
#'   \item Variable sequencing depth across spots
#'   \item Automatic calibration to a target mean count level
#'   \item Controlled zero proportion for non-SE genes
#' }
#'
#' @param spots A data frame or matrix of spatial coordinates
#'   with columns \code{x} and \code{y}.
#' @param num_genes Integer. Total number of genes to simulate.
#' @param prop Numeric vector of length 2 summing to 1.
#'   Proportions of non-SE and SE genes respectively
#'   (default: \code{c(0.8, 0.2)}).
#' @param pat_type Character vector specifying spatial pattern types for
#'   SE genes. Options include:
#'   \code{"hotspot"}, \code{"streak"}, \code{"gradient"},
#'   \code{"ring"}, \code{"wave"}, or \code{"mixed"}.
#'   If multiple values are supplied, patterns are assigned sequentially.
#' @param high_exp_prop Numeric. Proportion of spots in the
#'   high-expression region for spatial patterns (default: \code{0.3}).
#' @param grad_percent Numeric. Threshold controlling the spatial
#'   gradient region (default: \code{0.3}).
#' @param G Integer. Number of biological pathway groups.
#' @param gene_grp Integer vector of length \code{num_genes}
#'   specifying pathway membership for each gene.
#' @param rho Numeric. CAR spatial correlation parameter
#'   controlling dependence among genes within pathways
#'   (default: \code{0.85}).
#' @param tau_beta Numeric. CAR precision parameter governing
#'   variability of pathway effects (default: \code{10}).
#' @param eps_sd Numeric. Standard deviation of independent
#'   gene-level noise on the log scale (default: \code{0.20}).
#' @param eta_scale Numeric. Multiplicative scaling factor applied to
#'   spatial signal strength for SE genes (default: \code{3.0}).
#' @param max_count Integer. Maximum allowable expected count.
#'   Used to cap the baseline intensity during calibration
#'   (default: \code{500}).
#' @param zero_prop Numeric. Target proportion of zeros retained in
#'   non-SE genes after simulation. Excess zeros are replaced with
#'   small counts (default: \code{0.05}).
#' @param depth_model Character specifying the sequencing depth model.
#'   One of:
#'   \describe{
#'     \item{\code{"poisson"}}{Poisson-distributed sequencing depth}
#'     \item{\code{"negbin"}}{Negative binomial sequencing depth}
#'     \item{\code{"fixed"}}{Constant sequencing depth}
#'   }
#' @param depth_lambda Numeric. Mean sequencing depth parameter used for
#'   Poisson or negative binomial models (default: \code{5}).
#' @param depth_size Numeric. Dispersion parameter for the negative
#'   binomial depth model (default: \code{5}).
#' @param depth_fixed Integer. Fixed sequencing depth value when
#'   \code{depth_model = "fixed"} (default: \code{1}).
#' @param target_mean Numeric. Target median expected count used to
#'   automatically calibrate the baseline intensity parameter
#'   (default: \code{10}).
#' @param seed Integer or NULL. Random seed for reproducibility
#'   (default: NULL).
#' @param verbose Logical. If TRUE, prints a simulation summary
#'   including counts, depth statistics, and pattern distribution
#'   (default: TRUE).
#'
#' @return A named list containing:
#' \describe{
#'   \item{Y}{Integer matrix (\eqn{n \times p}) of simulated counts.}
#'   \item{spots}{Data frame of spatial coordinates used in simulation.}
#'   \item{N_i}{Integer vector of sequencing depth for each spot.}
#'   \item{mu0}{Auto-calibrated baseline log-mean expression level.}
#'   \item{Beta}{Numeric vector of CAR pathway effects.}
#'   \item{epsilon}{Matrix of independent noise terms.}
#'   \item{eta}{Matrix of spatial pattern effects for SE genes.}
#'   \item{G}{Number of biological pathway groups.}
#'   \item{Z}{Integer vector indicating gene state (1 = non-SE, 2 = SE).}
#'   \item{se_idx}{Integer vector of spatially expressed gene indices.}
#'   \item{gene_grp}{Integer vector of pathway membership.}
#'   \item{lambda}{Matrix of baseline Poisson intensities.}
#'   \item{Y_mean}{Matrix of expected counts before sampling.}
#'   \item{pat_types_used}{Character vector of assigned pattern types.}
#' }
#'
#' @details
#' Counts are generated according to a Poisson model:
#' \deqn{
#' Y_{ij} \sim \text{Poisson}(N_i \lambda_{ij}),
#' }
#' where sequencing depth \eqn{N_i} varies across spatial locations and
#' the log-intensity is defined as:
#' \deqn{
#' \log(\lambda_{ij}) =
#' \mu_0 + \beta_j + \eta_{ij} + \epsilon_{ij}.
#' }
#' The baseline parameter \eqn{\mu_0} is automatically calibrated so that
#' the simulated data achieve a target median count level specified by
#' \code{target_mean}.
#'
#' @importFrom stats rnorm rpois rnbinom median quantile
#'
#' @export
generate_genedata_pattern <- function(
    spots, num_genes, prop = c(0.8, 0.2),
    pat_type = "hotspot",
    high_exp_prop = 0.3, grad_percent = 0.3,
    G, gene_grp,
    rho = 0.85, tau_beta = 10,
    eps_sd = 0.20, eta_scale = 3.0,
    max_count = 500,
    zero_prop = 0.05,
    depth_model = c("poisson", "negbin", "fixed"),
    depth_lambda = 5, depth_size = 5, depth_fixed = 1,
    target_mean = 10,
    seed = NULL, verbose = TRUE
) {

  depth_model <- match.arg(depth_model)
  spots <- as.data.frame(spots)
  if (!all(c("x", "y") %in% names(spots))) names(spots)[1:2] <- c("x", "y")

  if (length(prop) != 2 || any(prop < 0) || abs(sum(prop) - 1) > 1e-8)
    stop("prop must be length-2 and sum to 1.")
  if (length(gene_grp) != num_genes)
    stop("gene_grp must have length num_genes.")

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(spots)
  P <- num_genes

  # --- Sequencing depth ---
  if (depth_model == "poisson") {
    N_i <- rpois(n, lambda = depth_lambda) + 1L
  } else if (depth_model == "negbin") {
    N_i <- rnbinom(n, mu = depth_lambda, size = depth_size) + 1L
  } else {
    N_i <- rep.int(as.integer(depth_fixed), n)
  }

  # --- Gene states ---
  n_se      <- round(P * prop[2])
  se_idx    <- sort(sample(seq_len(P), n_se))
  nonse_idx <- setdiff(seq_len(P), se_idx)
  Z <- rep(1L, P); Z[se_idx] <- 2L

  # --- CAR Beta ---
  Beta <- create_car_beta2(num_genes = P, gene_grp = gene_grp,
                           rho = rho, tau_beta = tau_beta)

  # --- Noise ---
  epsilon <- matrix(rnorm(n * P, 0, eps_sd), nrow = n, ncol = P)

  # --- Spatial eta for SE genes ---
  eta <- matrix(0, nrow = n, ncol = P)
  pat_types_used <- rep("NonSE", P)

  if (n_se > 0) {
    available <- c("hotspot", "streak", "gradient", "ring", "wave")
    if (length(pat_type) == 1 && pat_type == "mixed") {
      pat_vec <- rep_len(available, n_se)
    } else if (length(pat_type) == 1) {
      pat_vec <- rep(pat_type, n_se)
    } else {
      pat_vec <- rep_len(pat_type, n_se)
    }

    for (k in seq_along(se_idx)) {
      j <- se_idx[k]
      patt <- simulate_expression(spots, pat_vec[k], high_exp_prop, grad_percent, gene_id = j)
      eta[, j] <- patt$eta * eta_scale
      pat_types_used[j] <- pat_vec[k]
    }
  }

  # --- Auto-calibrate mu0 ---
  log_part <- sweep(epsilon + eta, 2, Beta, `+`)
  mu0      <- log(target_mean) - log(median(N_i)) - median(log_part)

  # --- Clamp mu0 so max expected count <= max_count ---
  mu0_cap <- log(max_count) - log(max(N_i)) - max(log_part)
  mu0     <- min(mu0, mu0_cap)

  # --- Generate counts ---
  lambda <- exp(mu0 + log_part)
  Y_mean <- outer(N_i, rep(1, P)) * lambda
  Y      <- matrix(rpois(n * P, lambda = as.vector(Y_mean)), nrow = n, ncol = P)

  # --- Replace excess zeros in Non-SE genes ---
  # Keep only zero_prop fraction as true zeros, rest become small counts (1-4)
  for (j in nonse_idx) {
    zero_idx <- which(Y[, j] == 0)
    n_zeros  <- length(zero_idx)
    if (n_zeros == 0) next

    # How many zeros to keep
    n_keep <- round(n * zero_prop)
    n_fill <- max(0, n_zeros - n_keep)

    if (n_fill > 0) {
      fill_idx <- sample(zero_idx, n_fill)
      Y[fill_idx, j] <- sample(1:4, n_fill, replace = TRUE)
    }
  }

  # --- Gene names ---
  gene_names <- character(P)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  for (j in se_idx) gene_names[j] <- paste0("Gene", j, "_SE_", pat_types_used[j])
  colnames(Y) <- gene_names

  if (verbose) {
    cat("---- ST Pattern Simulation Summary ----\n")
    cat("n spots:", n, " | P genes:", P, "\n")
    cat("SE genes:", n_se, " (", round(100 * n_se / P, 1), "%)\n", sep = "")
    if (n_se > 0) {
      pat_tab <- table(pat_types_used[se_idx])
      cat("Patterns:", paste(names(pat_tab), "=", pat_tab, collapse = ", "), "\n")
      cat("eta_scale:", eta_scale, "\n")
    }
    cat("zero_prop:", zero_prop, "\n")
    cat("Depth model:", depth_model,
        " | median N_i:", median(N_i),
        " | range:", paste(range(N_i), collapse = " - "), "\n")
    cat("mu0 (auto-calibrated):", signif(mu0, 4), "\n")
    cat("Non-SE zero %        :", round(100 * mean(Y[, nonse_idx] == 0), 1), "%\n")
    cat("Expected mean counts : median =", signif(median(Y_mean), 4),
        " | 95% =", paste(signif(quantile(as.vector(Y_mean), c(.025, .975)), 4), collapse = ", "),
        " | max =", signif(max(Y_mean), 4), "\n")
    cat("Observed counts      : median =", median(Y),
        " | 95% =", paste(quantile(as.vector(Y), c(.025, .975)), collapse = ", "),
        " | max =", max(Y), "\n")
    cat("---------------------------------------\n")
  }

  return(list(
    Y = Y, spots = spots, N_i = N_i, mu0 = mu0,
    Beta = Beta, epsilon = epsilon, eta = eta, G=G,
    Z = Z, se_idx = se_idx, gene_grp = gene_grp,
    lambda = lambda, Y_mean = Y_mean,
    pat_types_used = pat_types_used
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
