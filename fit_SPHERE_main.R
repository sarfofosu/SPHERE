## =============================================================================================

#' Fit the SPHERE Model for Spatial Transcriptomics Data
#'
#' Fits the SPHERE (Spatial Poisson Hierarchical modEl with pathway-infoRmed
#' gEne networks) model using a Bayesian framework implemented in Stan. The
#' model identifies spatially variable genes (SVGs) by combining a Poisson
#' log-normal likelihood, low-rank Gaussian Process spatial effects, and
#' pathway-informed Conditional Autoregressive (CAR) priors.
#'
#' @param data_mat A numeric matrix (\eqn{n \times p}) of gene expression
#'   counts, where \eqn{n} is the number of spatial spots and \eqn{p} is the
#'   number of genes. Non-integer values are rounded.
#' @param spot A numeric matrix (\eqn{n \times d}) of spatial coordinates for
#'   each spot (typically \eqn{d = 2} for x-y coordinates).
#' @param gene_group A character or factor vector of length \eqn{p} giving the
#'   pathway (group) membership for each gene.
#' @param stan_model_path Character string. Path to the compiled Stan model
#'   file (\code{.stan}).
#' @param iter_sampling Integer. Number of post-warmup sampling iterations per
#'   chain (default: 2000).
#' @param iter_warmup Integer. Number of warmup (burn-in) iterations per chain
#'   (default: 1000).
#' @param chains Integer. Number of MCMC chains (default: 3).
#' @param seed Integer. Random seed for reproducibility (default: 8).
#' @param knots Integer. Number of knots (inducing points) for the low-rank GP
#'   approximation (default: 30).
#' @param alpha Numeric vector of length 2. Dirichlet concentration parameters
#'   for the mixture weights (default: \code{c(10, 3)}). The first element
#'   corresponds to the non-spatially-variable component.
#' @param refresh Integer. Print progress every \code{refresh} iterations
#'   (default: 100). Set to 0 to suppress output.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{fit}{The \code{\link[cmdstanr]{CmdStanMCMC}} object containing the
#'     full posterior samples.}
#'   \item{summary}{A data frame of posterior summary statistics for key
#'     parameters (mean, median, sd, 95\% credible intervals, convergence
#'     diagnostics).}
#'   \item{draws}{A \code{draws_array} object from the \pkg{posterior} package
#'     containing all posterior draws.}
#'   \item{stan_data}{The list of data passed to Stan (useful for
#'     reproducibility and diagnostics).}
#'   \item{runtime}{Numeric. Total elapsed wall-clock time in seconds.}
#' }
#'
#' @details
#' The SPHERE model has three main components:
#' \enumerate{
#'   \item \strong{Observation model:} Counts follow a Poisson distribution
#'     with log-rates modeled as latent variables, offset by per-spot library
#'     size.
#'   \item \strong{Spatial structure:} A low-rank Radial Basis Function (RBF)
#'     Gaussian Process captures spatial expression patterns using \code{knots}
#'     inducing points.
#'   \item \strong{Gene classification:} A two-component mixture model with a
#'     gene-level indicator \eqn{Z_j} classifies each gene as spatially
#'     variable (SE) or non-spatially variable (non-SE). Pathway information
#'     is incorporated via a CAR prior on gene-level effects \eqn{\beta_j}.
#' }
#'
#' The posterior probability of spatial variability for gene \eqn{j} is
#' estimated as the fraction of posterior samples where \eqn{Z_j = 2}.
#'
#' @seealso \code{\link{make_rbf_basis}}, \code{\link{make_rbf_dist}}
#'
#' @examples
#' \dontrun{
#' result <- fit_sphere(
#'   data_mat = counts_matrix,
#'   spot = coord_matrix,
#'   gene_group = pathway_labels,
#'   stan_model_path = system.file("stan", "SPHERE.stan", package = "SPHERE"),
#'   knots = 50,
#'   chains = 4
#' )
#'
#' # Posterior probability of each gene being spatially variable
#' z_summary <- result$summary[grep("^Z\\[", result$summary$variable), ]
#' }
#'
#' @importFrom cmdstanr cmdstan_model
#' @importFrom posterior default_summary_measures default_convergence_measures
#'   quantile2
#' @importFrom stats dist median rnorm runif
#' @importFrom gtools rdirichlet
#' @export
#' 
#' 
## =============================================================================================

fit_sphere <- function(data_mat, spot, gene_group, stan_model_path, chains = 3, 
                       iter_sampling = 2000,  iter_warmup = 1000, knots = 30,
                       seed = 8, alpha = c(10, 3), refresh = 100) {
  
  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  if (!is.matrix(data_mat) && !is.data.frame(data_mat)) {
    stop("`data_mat` must be a matrix or data frame.", call. = FALSE)
  }
  data_mat <- round(as.matrix(data_mat))
  
  if (!is.matrix(spot) && !is.data.frame(spot)) {
    stop("`spot` must be a matrix or data frame of spatial coordinates.",
         call. = FALSE)
  }
  spot <- as.matrix(spot)
  
  if (nrow(data_mat) != nrow(spot)) {
    stop("Number of rows in `data_mat` (", nrow(data_mat),
         ") must match number of rows in `spot` (", nrow(spot), ").",
         call. = FALSE)
  }
  
  if (any(data_mat < 0)) {
    stop("`data_mat` must contain non-negative counts.", call. = FALSE)
  }
  
  # ------------------------------------------------------------------
  # 2. Extract dimensions
  # ------------------------------------------------------------------
  n <- nrow(data_mat)
  p <- ncol(data_mat)
  
  gene_grp <- as.integer(factor(gene_group))
  G <- length(unique(gene_grp))
  
  if (length(gene_grp) != p) {
    stop("`gene_group` must have length equal to the number of genes (",
         p, "), but has length ", length(gene_grp), ".", call. = FALSE)
  }
  
  if (knots >= n) {
    warning("Number of knots (", knots, ") >= number of spots (", n,
            "). Setting knots to n/2.", call. = FALSE)
    knots <- max(floor(n / 2), 1L)
  }
  
  # ------------------------------------------------------------------
  # 3. Normalization and spatial distances
  # ------------------------------------------------------------------
  N_i <- rowSums(data_mat)
  
  if (any(N_i == 0)) {
    stop("Some spots have zero total counts. Remove empty spots before ", "fitting.", call. = FALSE)
  }
  dist_sq <- as.matrix(dist(spot))^2
  
  # ------------------------------------------------------------------
  # 4. Low-rank GP basis construction
  # ------------------------------------------------------------------
  D   <- make_rbf_dist(spot, r = knots)
  Phi <- make_rbf_basis(spot, r = knots, lengthscale = NULL)
  r   <- knots
  
  # ------------------------------------------------------------------
  # 5. Assemble Stan data list
  # ------------------------------------------------------------------
  
  stan_data <- list(
    P = p,                                            # number of genes
    N = n,                                            # number of spatial locations
    Y = data_mat,                                     # observed count data
    dist_sq = dist_sq,                                # spatial distance matrix
    alpha = alpha,                                    # Hyperparameters: Dirichlet for mixture weights
    N_i = as.vector(N_i),                             # normalization/exposure term
    # Hyperparameters for priors
    a_err = 0, b_err = 1,                             # Hyperparameters: observation noise (half-normal)
    a_gsl = 0, b_gsl = 3,                             # Hyperparameters: GP lengthscale (log-normal)
    a_gs  = 0, b_gs  = 12,                            # Hyperparameters: GP amplitude (half-normal)
    a_rho = 2, b_rho = 2,                             # Hyperparameters: CAR correlation (beta)
    a_tau_beta = 1, b_tau_beta = 1,                   # Hyperparameters: CAR precision (half-normal)
    # Pathway structure
    G = G,                                            # number of pathways
    gene_group = gene_grp,                            # pathway membership per gene
    # Low-rank GP inputs
    r = r,
    D = D,
    Phi = Phi
  )
  
 
  # ------------------------------------------------------------------
  # 6. Initial values for MCMC (data-driven, robust to zeros)
  # ------------------------------------------------------------------
  
  # Compute safe log-rates from the data to initialize loglambda.
  # Replace zeros with 0.5 before taking log to avoid -Inf, then
  # subtract log(N_i) to get the rate on the correct scale.
  Y_safe      <- pmax(data_mat, 0.5)
  log_rates   <- log(Y_safe) - log(N_i)           # n x p matrix
  mu0_init    <- median(log_rates)                 # robust central value
  resid_init  <- log_rates - mu0_init              # gene-level deviations
  beta_init   <- colMeans(resid_init)              # per-gene mean residual
  ll_init     <- sweep(log_rates, 2, beta_init)    # loglambda ~ log_rate - beta
  med_dist    <- median(sqrt(D))                   # typical spot-knot distance
  
  init_fun <- function() {
    list(
      mu0        = mu0_init,
      pii        = replicate(p,as.numeric(gtools::rdirichlet(1, c(8, 2))), simplify = FALSE),
      loglambda  = ll_init,
      sigma_sd   = pmax(apply(log_rates, 2, sd), 0.1),
      sig_eta_gs = rep(0.3, p),
      ell_gs     = rep(med_dist, p),
      w          = matrix(rnorm(r * p, 0, 0.2), nrow = r, ncol = p),
      Beta       = beta_init,
      sigma_beta = 0.5,
      rho = runif(1, 0.7, 0.95)
    )
  }
  
  # ------------------------------------------------------------------
  # 7. Compile and fit model
  # ------------------------------------------------------------------

  # Compile Stan model from file
  model <- cmdstanr::cmdstan_model(stan_model_path)
  
  # Start timing model fitting
  t_start <- proc.time()[3]
  
  # Run MCMC sampling  
  fit <- model$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    seed            = seed,
    init            = init_fun,
    refresh         = refresh
  )
  
  # Compute total runtime  
  runtime <- proc.time()[3] - t_start
  
  message("SPHERE model fitting completed in ", round(runtime, 1), " seconds.")
  
  # ------------------------------------------------------------------
  # 8. Posterior summaries
  # ------------------------------------------------------------------
  
  # Extract summary statistics for key model parameters  
  summary <- fit$summary(
    variables = c("pii", "Z", "rho", "sig_eta_gs",
                  "Beta", "mu0", "sigma_beta", "ell_gs"),
    posterior::default_summary_measures()[1:3],
    quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)),
    posterior::default_convergence_measures()
  )
  
  # Extract posterior draws
  draws_array <- fit$draws()
  
  # ------------------------------------------------------------------
  # 9. Return results
  # ------------------------------------------------------------------
  structure(
    list(
      fit = fit,                                     # full CmdStan object
      summary = summary,                             # posterior summaries
      draws = draws_array,                           # selected posterior samples
      stan_data = stan_data,                         # stan data 
      runtime = runtime                              # computation time    
    ),
    class = "sphere_fit"
  )
}
































