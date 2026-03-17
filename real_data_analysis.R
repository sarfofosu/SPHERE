
na_idx <- which(is.na(pathway_df$Pathway))       # Identify genes with missing pathway annotations

# If missing pathways exist, assign them into 3 artificial groups (p1, p2, p3)
if (length(na_idx) > 0) {     
  groups <- split(            # This prevents loss of genes and allows them to be included in CAR structure
    na_idx,
    cut(seq_along(na_idx), 3, labels = c("p1","p2","p3"))
  )
  
  # Replace NA pathway labels with artificial group labels
  for (g in names(groups)) {
    pathway_df$Pathway[groups[[g]]] <- g
  }
}





# Filter pathways based on size/criteria (user-defined helper function)
gen_path <- filter_pathways_by_limit(pathway_df)




fit_sphere <- function(data_mat, spot, pathway_df, stan_model_path, 
                       iter_sampling = 5000, iter_warmup = 2000,
                       chains = 3, seed = 8) {
  
    # ----------------------------
    # 1. Data preprocessing
    # ----------------------------
    data_mat <- round(as.matrix(data_mat))           # Ensure input data is a numeric matrix of counts (Poisson requires integers)
    na_idx <- which(is.na(pathway_df$Pathway))       # Identify genes with missing pathway annotations

    # If missing pathways exist, assign them into 3 artificial groups (p1, p2, p3)
    if (length(na_idx) > 0) {     
      groups <- split(            # This prevents loss of genes and allows them to be included in CAR structure
        na_idx,
        cut(seq_along(na_idx), 3, labels = c("p1","p2","p3"))
      )
      
      # Replace NA pathway labels with artificial group labels
      for (g in names(groups)) {
        pathway_df$Pathway[groups[[g]]] <- g
      }
    }
    
    # Filter pathways based on size/criteria (user-defined helper function)
    gen_path <- filter_pathways_by_limit(pathway_df)
    
    # ----------------------------
    # 2. Extract dimensions
    # ----------------------------
    n <- nrow(data_mat)                               # n = number of spatial locations (spots)
    p <- ncol(data_mat)                               # p = number of genes
    G <- length(unique(gen_path$Pathway))             # G = number of unique biological pathways
    gene_grp <- as.integer(factor(gen_path$Pathway))  # Convert pathway labels into integer indices for Stan
    
    # ----------------------------
    # 3. Construct model inputs
    # ----------------------------
    
    # Compute normalization factor E_ij:
    # total counts per spot replicated across genes
    # (acts like exposure/offset in Poisson model)
    E_ij <- matrix(rep(apply(data_mat, 1, sum), p), n, p)
    
    # Compute squared Euclidean distance matrix between spatial locations
    # Used in Gaussian process covariance
    dist_sq <- as.matrix(dist(spot))^2
    
    # ----------------------------
    # 4. Low-rank GP basis construction (NEW)
    # ----------------------------
    
    # Construct distance matrix from spots to RBF knots
    # D[i,b] = squared distance from spot i to knot b
    D <- make_rbf_dist(spot, r = 30)
    
    # Number of basis functions (knots)
    r <- ncol(D)
    
    # Median distance (can be useful for diagnostics / scaling)
    med_D <- median(D)
    
    # Construct RBF basis matrix
    # Phi[i,b] = basis function linking spot i to knot b
    Phi <- make_rbf_basis(spot, r = 30, lengthscale = NULL)
    
    # Ensure consistency
    r <- ncol(Phi)
    
    
    
    
    # Assemble data list to pass into Stan model
    stan_data <- list(
      P = p,                                            # number of genes
      N = n,                                            # number of spatial locations
      Y = data_mat,                                     # observed count data
      dist_sq = dist_sq,                                # spatial distance matrix
      alpha = c(10, 3),                                 # Dirichlet prior for mixture weights
      N_i = as.vector(N_i),                                     # normalization/exposure term
      
      # Hyperparameters for priors
      a_err = 0, b_err = 1,                             # half-normal prior for noise sd
      a_gsl = 0, b_gsl = 3,                             # lognormal prior for GP length-scale
      a_gs  = 0, b_gs  = 12,                            # half-normal prior for GP variance
      a_rho = 2, b_rho = 2,                             # beta prior for CAR correlation
      a_tau_beta = 1, b_tau_beta = 1,                   # gamma prior for CAR precision
      
      # Pathway structure
      G = G,                                            # number of pathways
      gene_group = gene_grp,                             # pathway membership per gene
      
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
        w = matrix(rnorm(r * p), r, p),                 # NEW: GP weights (r x p matrix)
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



