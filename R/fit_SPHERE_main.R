## =============================================================================================

#' Fit the SPHERE Model for Spatial Transcriptomics Data
#'
#' Fits the SPHERE (Spatial Poisson Hierarchical modEl with pathway-infoRmed
#' gEne networks) model using a Bayesian framework implemented in Stan. The
#' model identifies spatially expressed (SE) genes by combining a Poisson
#' log-normal likelihood, a low-rank Gaussian Process spatial effect, and a
#' pathway-informed Conditional Autoregressive (CAR) prior on each gene's SE
#' log-odds.
#'
#' @param data_mat A numeric matrix (\eqn{n \times p}) of gene expression
#'   counts, where \eqn{n} is the number of spatial spots and \eqn{p} is the
#'   number of genes. Non-integer values are rounded.
#' @param spot A numeric matrix or data frame (\eqn{n \times 2}) of spatial
#'   coordinates for each spot.
#' @param gene_group A vector of length \eqn{p} giving a single pathway label
#'   for each gene. Supply this or \code{pathway_membership}, not both
#'   (default: NULL).
#' @param pathway_membership A binary matrix (\eqn{p \times K}) of gene-pathway
#'   membership, allowing genes to belong to several pathways. Supply this or
#'   \code{gene_group}, not both (default: NULL).
#' @param N_i A numeric vector of length \eqn{n} giving the library size
#'   (sequencing depth) for each spot. If NULL, defaults to the row sums
#'   of \code{data_mat} (default: NULL).
#' @param inits A function returning a named list of initial values for the
#'   MCMC sampler, or NULL to use the default data-driven initialization
#'   (default: NULL).
#' @param stan_model_path Character string. Path to the Stan model file
#'   (\code{.stan}). Defaults to the model bundled with the package.
#' @param use_pathway Logical. If TRUE, the pathway network informs SE
#'   detection; if FALSE, every gene shares the same prior SE log-odds
#'   (default: TRUE).
#' @param fix_rho_zero Logical. If TRUE, the CAR dependence \eqn{\rho_u} is
#'   fixed at 0; if FALSE, it is estimated (default: FALSE).
#' @param car_type Character. Pathway CAR prior: \code{"leroux"} (Leroux CAR)
#'   or \code{"sparse"} (sparse proper CAR) (default: \code{"leroux"}).
#' @param use_shared_Phi Logical. If TRUE, one RBF basis is shared by all genes
#'   and the GP basis weights are integrated out analytically; if FALSE, each
#'   gene has its own lengthscale and the basis weights are sampled
#'   (default: TRUE).
#' @param iter_sampling Integer. Number of post-warmup sampling iterations per
#'   chain (default: 2000).
#' @param iter_warmup Integer. Number of warmup (burn-in) iterations per chain
#'   (default: 1000).
#' @param chains Integer. Number of MCMC chains (default: 3).
#' @param seed Integer. Random seed for knot selection and sampling
#'   (default: 8).
#' @param knots Integer. Number of knots (inducing points) for the low-rank GP
#'   approximation (default: 30).
#' @param refresh Integer. Print progress every \code{refresh} iterations
#'   (default: 100). Set to 0 to suppress output.
#' @param keep_fit Logical. If TRUE, the full \code{CmdStanMCMC} object is
#'   returned in \code{fit}. Keep FALSE for large numbers of genes, since the
#'   fit caches every variable that is loaded from it (default: FALSE).
#' @param save_draws Character vector. Variables whose posterior draws are
#'   returned in \code{draws}. Variables not present in the fitted model (e.g.
#'   \code{rho_u} when \code{fix_rho_zero = TRUE}) are skipped
#'   (default: \code{c("rho_u", "sigma_u", "alpha0", "mu0", "sigma_beta",
#'   "prob_SE", "sig_eta_gs")}).
#' @param mu_noise Numeric. Location of the half-normal prior on the
#'   gene-specific residual SD (default: 0).
#' @param sd_noise Numeric. Scale of the half-normal prior on the
#'   gene-specific residual SD (default: 0.5).
#' @param mu_gp_amplitude Numeric. Location of the half-normal prior on the GP
#'   amplitude (default: 0).
#' @param sd_gp_amplitude Numeric. Scale of the half-normal prior on the GP
#'   amplitude (default: 0.4).
#' @param mu_gp_lengthscale Numeric or NULL. Log-scale mean of the log-normal
#'   prior on the GP lengthscale. If NULL, \code{log(0.25 * domain_scale)},
#'   where \code{domain_scale} is the square root of the area spanned by
#'   \code{spot} (default: NULL).
#' @param sd_gp_lengthscale Numeric. Log-scale SD of the log-normal prior on
#'   the GP lengthscale (default: 0.4).
#' @param mu_intercept Numeric or NULL. Mean of the normal prior on the global
#'   intercept \eqn{\mu_0}. If NULL, the median log rate of the data
#'   (default: NULL).
#' @param sd_intercept Numeric. SD of the normal prior on \eqn{\mu_0}
#'   (default: 1).
#' @param mu_beta_sig Numeric. Location of the half-normal prior on the SD of
#'   the gene-level effects \eqn{\beta_j} (default: 0).
#' @param sd_beta_sig Numeric. Scale of the half-normal prior on the SD of the
#'   gene-level effects \eqn{\beta_j} (default: 0.5).
#' @param alpha0_mean Numeric. Mean of the normal prior on the global SE
#'   log-odds \eqn{\alpha_0} (default: \code{qlogis(0.2)}).
#' @param alpha0_sd Numeric. SD of the normal prior on \eqn{\alpha_0}
#'   (default: 1).
#' @param mu_u_sig Numeric. Location of the half-normal prior on the SD of the
#'   pathway effects \eqn{u_j} (default: 0).
#' @param sd_u_sig Numeric. Scale of the half-normal prior on the SD of the
#'   pathway effects \eqn{u_j} (default: 1.5).
#' @param shape_rho_u Numeric. First shape parameter of the Beta prior on the
#'   CAR dependence \eqn{\rho_u} (default: 1).
#' @param rate_rho_u Numeric. Second shape parameter of the Beta prior on the
#'   CAR dependence \eqn{\rho_u} (default: 1).
#'
#' @return A named list of class \code{"sphere_fit"} with the following
#'   elements:
#' \describe{
#'   \item{fit}{The \code{CmdStanMCMC} object containing the
#'     full posterior samples, or NULL when \code{keep_fit = FALSE}.}
#'   \item{summary}{A data frame of posterior summary statistics for key
#'     parameters (mean, median, sd, 95\% credible intervals, convergence
#'     diagnostics), with gene names on gene-level rows.}
#'   \item{prob_SE}{Named numeric vector of posterior mean probabilities that
#'     each gene is spatially expressed.}
#'   \item{draws}{A \code{draws_array} object from the \pkg{posterior} package
#'     containing the posterior draws of the variables in \code{save_draws}
#'     that exist in the fitted model.}
#'   \item{sampler_diag}{Sampler diagnostics (divergences, maximum treedepth
#'     hits, E-BFMI) from \code{fit$diagnostic_summary()}.}
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
#'   \item \strong{Gene classification:} A gene-level indicator \eqn{Z_j}
#'     classifies each gene as SE or non-SE. Its prior probability is
#'     \eqn{\mathrm{logit}^{-1}(\alpha_0 + u_j)}, where the pathway effects
#'     \eqn{u_j} follow a CAR prior on the pathway network: genes that share a
#'     pathway are neighbours.
#' }
#'
#' The posterior probability that gene \eqn{j} is SE is returned in
#' \code{prob_SE}. The same call fits real data and data simulated with either
#' \code{\link{generate_genedata_model}} (model-based) or
#' \code{\link{make_pattern_design}} (pattern-based); the input data alone
#' determine what is fitted.
#'
#' Model variants used in the paper:
#' \itemize{
#'   \item \code{use_pathway = TRUE, fix_rho_zero = FALSE}: pathway-informed
#'     model (M1).
#'   \item \code{use_pathway = TRUE, fix_rho_zero = TRUE}: \eqn{\rho_u = 0}
#'     (M2).
#'   \item \code{use_pathway = FALSE}: no pathway information (M3).
#' }
#'
#' If \code{N_i} is provided, it overrides the default row-sum normalization.
#' This is useful when library sizes have been pre-computed or adjusted
#' externally (e.g. after normalization pipelines).
#'
#' If \code{inits} is provided, it must be a function with no arguments that
#' returns a named list matching the Stan parameter block. If NULL, a
#' data-driven initialization is used automatically.
#'
#' @seealso \code{\link{generate_genedata_model}},
#'   \code{\link{make_sim_design}}, \code{\link{make_pattern_design}},
#'   \code{\link{simulate_counts}}
#'
#' @examples
#' \dontrun{
#' # Model-based simulated data
#' M   <- generate_pathway_membership(num_genes = 50, num_pathways = 10, seed = 1)
#' sim <- generate_genedata_model(num_spots = 100, num_genes = 50,
#'                                pathway_membership = M, seed = 2)
#' result <- fit_sphere(
#'   data_mat           = sim$Y,
#'   spot               = sim$spots,
#'   pathway_membership = sim$pathway_membership,
#'   N_i                = sim$N_i
#' )
#'
#' # Pattern-based simulated data: same call
#' design <- make_pattern_design(num_spots = 100, num_genes = 50,
#'                               pathway_membership = M, design_seed = 3)
#' pat    <- simulate_counts(design, seed = 4)
#' result <- fit_sphere(
#'   data_mat           = pat$Y,
#'   spot               = pat$spots,
#'   pathway_membership = pat$pathway_membership,
#'   N_i                = pat$N_i
#' )
#'
#' # Real data with one pathway label per gene, sparse CAR
#' result <- fit_sphere(
#'   data_mat   = counts_matrix,
#'   spot       = coord_matrix,
#'   gene_group = pathway_labels,
#'   car_type   = "sparse"
#' )
#'
#' # Posterior probability of each gene being spatially expressed
#' result$prob_SE
#' }
#'
#' @importFrom cmdstanr cmdstan_model
#' @importFrom posterior default_summary_measures default_convergence_measures
#'   quantile2
#' @importFrom stats median rnorm runif sd qlogis
#' @export

## =============================================================================================

fit_sphere <- function(data_mat, spot, gene_group = NULL, pathway_membership = NULL,
                       N_i = NULL, inits = NULL,
                       stan_model_path = system.file("stan", "SPHERE_stan.stan", package = "SPHERE"),
                       use_pathway = TRUE, fix_rho_zero = FALSE,
                       car_type = c("leroux", "sparse"), use_shared_Phi = TRUE,
                       iter_sampling = 2000, iter_warmup = 1000,
                       chains = 3, seed = 8, knots = 30, refresh = 100,
                       keep_fit = FALSE,
                       save_draws = c("rho_u", "sigma_u", "alpha0", "mu0",
                                      "sigma_beta", "prob_SE", "sig_eta_gs"),
                       mu_noise = 0, sd_noise = 0.5,
                       mu_gp_amplitude = 0, sd_gp_amplitude = 0.4,
                       mu_gp_lengthscale = NULL, sd_gp_lengthscale = 0.4,
                       mu_intercept = NULL, sd_intercept = 1,
                       mu_beta_sig = 0, sd_beta_sig = 0.5,
                       alpha0_mean = qlogis(0.2), alpha0_sd = 1,
                       mu_u_sig = 0, sd_u_sig = 1.5,
                       shape_rho_u = 1, rate_rho_u = 1) {

  car_type <- match.arg(car_type)

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
  spot <- as.matrix(spot)[, 1:2, drop = FALSE]

  if (nrow(data_mat) != nrow(spot)) {
    stop("Number of rows in `data_mat` (", nrow(data_mat),
         ") must match number of rows in `spot` (", nrow(spot), ").",
         call. = FALSE)
  }

  if (any(data_mat < 0)) {
    stop("`data_mat` must contain non-negative counts.", call. = FALSE)
  }

  if (!is.null(inits) && !is.function(inits)) {
    stop("`inits` must be a function or NULL.", call. = FALSE)
  }

  if (use_pathway && is.null(gene_group) && is.null(pathway_membership)) {
    stop("Supply `gene_group` or `pathway_membership`, or set use_pathway = FALSE.",
         call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 2. Extract dimensions
  # ------------------------------------------------------------------
  n <- nrow(data_mat)
  p <- ncol(data_mat)

  gene_names <- colnames(data_mat)
  if (is.null(gene_names)) gene_names <- paste0("Gene", seq_len(p))

  if (knots >= n) {
    warning("Number of knots (", knots, ") >= number of spots (", n,
            "). Setting knots to n/2.", call. = FALSE)
    knots <- max(floor(n / 2), 1L)
  }
  r <- knots

  if (!use_pathway) fix_rho_zero <- TRUE                 # no pathway: rho_u is not used

  # ------------------------------------------------------------------
  # 3. Normalization (library sizes)
  # ------------------------------------------------------------------
  if (is.null(N_i)) {
    N_i <- rowSums(data_mat)
    message("N_i not provided. Using row sums of data_mat as library sizes.")
  } else {
    N_i <- as.numeric(N_i)
    if (length(N_i) != n) {
      stop("`N_i` must have length equal to the number of spots (",
           n, "), but has length ", length(N_i), ".", call. = FALSE)
    }
    message("Using user-supplied N_i as library sizes.")
  }

  if (any(N_i <= 0)) {
    stop("Some spots have zero total counts. Remove empty spots before ",
         "fitting.", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 4. Pathway network for the CAR prior
  # ------------------------------------------------------------------
  # Genes that share at least one pathway are neighbours.
  # Without pathway information, an empty network is passed to Stan.
  if (!is.null(gene_group) || !is.null(pathway_membership)) {
    W <- build_pathway_network(p, gene_grp = gene_group,
                               pathway_membership = pathway_membership)$W
  } else {
    W <- matrix(0L, p, p)
  }
  car_data <- prepare_sparse_car_W(W)

  # ------------------------------------------------------------------
  # 5. Low-rank GP basis construction
  # ------------------------------------------------------------------
  # Knots, distances and basis come from one seeded k-means, so D and Phi match.
  basis <- make_rbf_knots(spot, r = r, seed = seed)
  D     <- basis$D
  Phi   <- basis$Phi

  # ------------------------------------------------------------------
  # 6. Data-driven defaults and initial values
  # ------------------------------------------------------------------
  Y_safe        <- pmax(data_mat, 0.5)
  log_rates     <- log(Y_safe) - log(N_i)                # n x p matrix of log-rates
  mu0_init      <- median(log_rates)                     # robust global intercept
  resid_init    <- log_rates - mu0_init                  # gene-level deviations
  beta_init     <- pmax(pmin(colMeans(resid_init), 1), -1)
  ll_init       <- log_rates                             # initial loglambda
  med_dist      <- median(D[D > 0])                      # typical spot-knot distance
  ell_init      <- 0.5 * med_dist
  sigma_sd_init <- pmax(apply(log_rates, 2, sd), 0.1)

  domain_scale <- sqrt(diff(range(spot[, 1])) * diff(range(spot[, 2])))
  if (is.null(mu_gp_lengthscale)) mu_gp_lengthscale <- log(0.25 * domain_scale)
  if (is.null(mu_intercept))      mu_intercept      <- mu0_init

  if (is.null(inits)) {
    # Default data-driven initialization
    init_fun <- function() {
      mu0_i      <- mu0_init + rnorm(1, 0, 0.2)
      Beta_i     <- beta_init + rnorm(p, 0, 0.05)
      sigma_sd_i <- sigma_sd_init * runif(1, 0.8, 1.25)
      init <- list(
        mu0        = mu0_i,
        z_ll       = sweep(sweep(ll_init, 2, mu0_i + Beta_i), 2, sigma_sd_i, "/"),
        sigma_sd   = sigma_sd_i,
        sig_eta_gs = rep(0.20, p) * runif(1, 0.7, 1.4),
        ell_gs     = rep(ell_init, p) * runif(p, 0.8, 1.25),
        Beta       = Beta_i,
        sigma_beta = runif(1, 0.2, 0.6),
        alpha0     = alpha0_mean + rnorm(1, 0, 0.1)
      )
      # Zero-sized Stan parameters are not initialized.
      if (!use_shared_Phi) {
        init$w <- matrix(rnorm(r * p, 0, 0.2), nrow = r, ncol = p)
      }
      if (use_pathway) {
        init$u       <- rnorm(p, 0, 0.1)
        init$sigma_u <- array(runif(1, 0.5, 1.5), dim = 1)
      }
      if (use_pathway && !fix_rho_zero) {
        init$rho_u_free <- array(runif(1, 0.35, 0.8), dim = 1)
      }
      init
    }
  } else {
    # Use user-supplied initialization function
    init_fun <- inits
  }

  # ------------------------------------------------------------------
  # 7. Assemble Stan data list
  # ------------------------------------------------------------------
  stan_data <- list(
    # Data dimensions and observed data
    N   = n,                                             # number of spatial spots
    P   = p,                                             # number of genes
    Y   = data_mat,                                      # observed count matrix (n x p)
    N_i = as.vector(N_i),                                # per-spot library size (exposure)
    # Low-rank GP inputs
    r   = r,                                             # number of knots
    Phi = Phi,                                           # shared RBF basis (n x r)
    D   = D,                                             # spot-to-knot distances (n x r)
    # Model options
    use_shared_Phi = as.integer(use_shared_Phi),
    use_pathway    = as.integer(use_pathway),
    fix_rho_zero   = as.integer(fix_rho_zero),
    use_car_lx     = as.integer(car_type == "leroux"),   # 1 = Leroux CAR; 0 = sparse proper CAR
    K_u            = if (use_pathway) p else 0L,
    K_sigma_u      = if (use_pathway) 1L else 0L,
    K_rho          = if (use_pathway && !fix_rho_zero) 1L else 0L,
    # Residual noise prior (half-normal)
    mu_noise = mu_noise, sd_noise = sd_noise,
    # GP amplitude prior (half-normal)
    mu_gp_amplitude = mu_gp_amplitude, sd_gp_amplitude = sd_gp_amplitude,
    # GP lengthscale prior (log-normal)
    mu_gp_lengthscale = mu_gp_lengthscale, sd_gp_lengthscale = sd_gp_lengthscale,
    # Global intercept prior (normal)
    mu_intercept = mu_intercept, sd_intercept = sd_intercept,
    # Pathway / CAR structure
    P_car     = car_data$P_car,
    car_index = car_data$car_index,
    P_iso     = car_data$P_iso,
    iso_index = car_data$iso_index,
    W_n       = car_data$W_n,
    W_sparse  = car_data$W_sparse,
    D_sparse  = car_data$D_sparse,
    lap_eig   = car_data$lap_eig,                        # Leroux CAR log determinant
    car_eig   = car_data$car_eig,                        # sparse CAR log determinant
    # Gene-level effect prior (half-normal on sigma_beta)
    mu_beta_sig = mu_beta_sig, sd_beta_sig = sd_beta_sig,
    # Pathway-informed detection priors
    alpha0_mean = alpha0_mean, alpha0_sd = alpha0_sd,
    mu_u_sig = mu_u_sig, sd_u_sig = sd_u_sig,
    shape_rho_u = shape_rho_u, rate_rho_u = rate_rho_u
  )

  # ------------------------------------------------------------------
  # 8. Compile and fit model
  # ------------------------------------------------------------------
  if (!file.exists(stan_model_path)) {
    stop("Stan model file not found at: ", stan_model_path,
         "\nRe-install the SPHERE package.", call. = FALSE)
  }

  model   <- cmdstanr::cmdstan_model(stan_model_path)
  set.seed(seed)
  t_start <- proc.time()[3]

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

  runtime <- proc.time()[3] - t_start
  message("SPHERE model fitting completed in ", round(runtime, 1), " seconds.")

  # ------------------------------------------------------------------
  # 9. Posterior summaries
  # ------------------------------------------------------------------
  summary_vars <- c("alpha0", "pi", "Z", "prob_SE", "logLR", "sig_eta_gs",
                    "Beta", "sigma_sd", "mu0", "sigma_beta")
  if (use_pathway) summary_vars <- c("sigma_u", "u", "u_centered", summary_vars)
  if (use_pathway && !fix_rho_zero) summary_vars <- c("rho_u", summary_vars)
  if (!use_shared_Phi) summary_vars <- c(summary_vars, "ell_gs")

  summary <- fit$summary(
    variables = summary_vars,
    posterior::default_summary_measures()[1:3],
    quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)),
    posterior::default_convergence_measures()
  )

  # Add gene names to gene-level rows
  for (v in c("Z", "prob_SE", "logLR", "Beta", "sig_eta_gs", "sigma_sd",
              "u", "u_centered", "pi", "ell_gs")) {
    idx <- grep(paste0("^", v, "\\["), summary$variable)
    if (length(idx) == p) summary$variable[idx] <- paste0(v, "[", gene_names, "]")
  }

  # Posterior probability that each gene is SE
  prob_SE <- as.numeric(colMeans(fit$draws(variables = "prob_SE", format = "draws_matrix")))
  names(prob_SE) <- gene_names

  # Keep draws only for the requested variables that exist in this model
  draw_vars    <- intersect(save_draws, fit$metadata()$stan_variables)
  draws_array  <- if (length(draw_vars) > 0) fit$draws(variables = draw_vars) else NULL
  sampler_diag <- fit$diagnostic_summary()

  # ------------------------------------------------------------------
  # 10. Return results
  # ------------------------------------------------------------------
  structure(
    list(
      fit          = if (keep_fit) fit else NULL,  # full CmdStanMCMC object (optional)
      summary      = summary,       # posterior summary statistics
      prob_SE      = prob_SE,       # posterior P(gene is SE)
      draws        = draws_array,   # posterior draws of save_draws variables
      sampler_diag = sampler_diag,  # divergences, treedepth, E-BFMI
      stan_data    = stan_data,     # data passed to Stan
      runtime      = runtime        # total wall-clock time in seconds
    ),
    class = "sphere_fit"
  )
}
