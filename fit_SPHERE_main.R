#' Fit the SPHERE Model for Spatial Transcriptomics Data
#'
#' This function fits the SPHERE (A Spatial Poisson Hierarchical modEl with pathway-infoRmed gEne networks) 
#' using a Bayesian framework implemented in Stan.
#'
#' @param data_mat A numeric matrix (n x p) of gene expression counts,
#'   where n is the number of spatial locations (spots) and p is the number of genes.
#' @param spot A numeric matrix (n x d) of spatial coordinates for each spot.
#' @param pathway_df A data frame containing gene-to-pathway mappings.
#'   Must include a column named `Pathway`.
#' @param stan_model_path Path to the Stan model file.
#' @param iter_sampling Number of post-warmup iterations (default: 2000).
#' @param iter_warmup Number of warmup iterations (default: 1000).
#' @param chains Number of MCMC chains (default: 1).
#' @param seed Random seed for reproducibility (default: 8).
#'
#' @return A list containing:
#' \describe{
#'   \item{fit}{CmdStanMCMC object}
#'   \item{summary}{Posterior summary statistics}
#'   \item{runtime}{Elapsed time (seconds)}
#'   \item{draws}{Posterior draws (selected parameters)}
#' }
#'
#' @details
#' The model:
#' \itemize{
#'   \item Uses a Poisson likelihood for count data
#'   \item Incorporates spatial dependence via Gaussian processes
#'   \item Models gene relationships using pathway-informed CAR priors
#' }
#'
#' @examples
#' \dontrun{
#' result <- fit_sphere(
#'   data_mat = counts,
#'   spot = coordinates,
#'   pathway_df = pathway_data,
#'   stan_model_path = "model_w_fullCAR.stan"
#' )
#' }
#'
#' @export
fit_sphere <- function(data_mat, spot, pathway_df, stan_model_path, 
                       iter_sampling = 5000, iter_warmup = 2000,
                       chains = 3, gene_group, seed = 8, knots = 30) {
  
  # ----------------------------
  # 1. Data
  # ----------------------------
  data_mat <- round(as.matrix(data_mat))           # Ensure input data is a numeric matrix of counts (Poisson requires integers)

  # ----------------------------
  # 2. Extract dimensions
  # ----------------------------
  n <- nrow(data_mat)                               # n = number of spatial locations (spots)
  p <- ncol(data_mat)                               # p = number of genes
  G <- length(unique(gene_group))                   # G = number of unique biological pathways
  gene_grp <- as.integer(factor(gene_group))        # Convert pathway labels into integer indices for Stan
  
  if (length(gene_grp) != p) {
    stop("gene_group must have length equal to number of genes (columns of data_mat).")
  }
  
  if (max(gene_grp) > G) {
    stop("gene_group contains values larger than G.")
  }
  
  # ----------------------------
  # 3. Construct model inputs
  # ----------------------------
  
  # Compute normalization factor N_i: total counts per spot replicated across genes
  N_i <- apply(data_mat, 1, sum)
  # Compute squared Euclidean distance matrix between spatial locations
  dist_sq <- as.matrix(dist(spot))^2
  
  # ----------------------------
  # 4. Low-rank GP basis construction
  # ----------------------------
  
  # Construct distance matrix from spots to RBF knots
  D <- make_rbf_dist(spot, r = knots)           # D[i,b] = squared distance from spot i to knot b
  med_D <- median(D)                            # Median distance (can be useful for diagnostics / scaling)
  # Construct RBF basis matrix
  Phi <- make_rbf_basis(spot, r = knots, lengthscale = NULL)   # Phi[i,b] = basis function linking spot i to knot b
  r <- knots                                    # Number of basis functions (knots)
  
  # Assemble data list to pass into Stan model
  stan_data <- list(
    P = p,                                            # number of genes
    N = n,                                            # number of spatial locations
    Y = data_mat,                                     # observed count data
    dist_sq = dist_sq,                                # spatial distance matrix
    alpha = c(10, 3),                                 # Dirichlet prior for mixture weights
    N_i = as.vector(N_i),                             # normalization/exposure term
    
    # Hyperparameters for priors
    a_err = 0, b_err = 1,                             # half-normal prior for noise sd
    a_gsl = 0, b_gsl = 3,                             # lognormal prior for GP length-scale
    a_gs  = 0, b_gs  = 12,                            # half-normal prior for GP variance
    a_rho = 2, b_rho = 2,                             # beta prior for CAR correlation
    a_tau_beta = 1, b_tau_beta = 1,                   # gamma prior for CAR precision
    
    # Pathway structure
    G = G,                                            # number of pathways
    gene_group = gene_grp,                            # pathway membership per gene
    
    # Low-rank GP inputs
    r = r,
    D = D,
    Phi = Phi
  )
  
  # ----------------------------
  # 4. Initial values for MCMC
  # ----------------------------
  init_fun <- function() {
    list(
      sig_eta_gs = rep(10, p),                        # Gene-specific GP variance parameters
      ell_gs = rep(2, p),                             # Gene-specific length-scale parameters
      pii = replicate(                                # Mixture probabilities (Dirichlet initialized)
        p, as.numeric(gtools::rdirichlet(1, c(8, 2))), simplify = FALSE),
      sigma_sd = rep(1, p),                           # Observation noise standard deviation
      Beta = rep(1, p),                               # CAR regression coefficients
      rho = runif(1, 0.1, 1),                         # Spatial correlation parameter
      sigma_beta = runif(1, 0.1, 1),                  # CAR precision parameter
      mu0 = 1,                                        # Global intercept
      loglambda = matrix(0.5, n, p),                  # Log-intensity (Poisson mean parameter)
      w = matrix(rnorm(r * p), r, p),                 # GP weights (r x p matrix)
    )
  }
  
  # ----------------------------
  # 5. Compile and fit model
  # ----------------------------
  
  # Compile Stan model from file
  model <- cmdstanr::cmdstan_model(stan_model_path)
  
  # Start timing model fitting
  t_start <- proc.time()[3]
  
  # Run MCMC sampling
  fit <- model$sample(
    data = stan_data,                                # input data
    chains = chains,                                 # number of chains
    parallel_chains = chains,                        # parallel execution
    iter_warmup = iter_warmup,                       # burn-in iterations
    iter_sampling = iter_sampling,                   # posterior samples
    seed = seed,                                     # reproducibility
    init = init_fun,                                 # initialization function
    refresh = 100                                    # print progress every 100 iterations
  )
  
  # Compute total runtime
  runtime <- proc.time()[3] - t_start
  
  # ----------------------------
  # 6. Posterior summaries
  # ----------------------------
  
  # Extract summary statistics for key model parameters
  summary <- fit$summary(
    variables = c("pii", "Z", "rho", "sig_eta_gs", "Beta", "mu0","sigma_beta", "ell_gs"),
    posterior::default_summary_measures()[1:3],
    quantiles = ~ quantile2(., probs = c(0.025, 0.975)),
    posterior::default_convergence_measures()
  )    
  
  
  ##-----------------------------------------------
  ## Extract posterior draws
  ##-----------------------------------------------
  
  draws_array <- fit$draws()
  
  # ----------------------------
  # 8. Return results
  # ----------------------------
  
  return(list(
    fit = fit,                                     # full CmdStan object
    summary = summary,                             # posterior summaries
    draws = draws_array,                           # selected posterior samples
    stan_data = stan_data,                         # stan data 
    runtime = runtime))                            # computation time    
}



