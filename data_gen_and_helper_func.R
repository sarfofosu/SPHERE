

## ------------------------------------------------------------
## Function to generate fixed number of spots (coordinates)
## ------------------------------------------------------------

generate_spatial_spots <- function(num_spots, x1 = 0, x2 = 1, y1 = 0, y2 = 1, seed =NULL) {
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data.frame(
    x = runif(num_spots, min = x1, max = x2),
    y = runif(num_spots, min = y1, max = y2)
  )
}


## ------------------------------------------------------------
## Function to simulate a name spatial pattern
## ------------------------------------------------------------

simulate_expression <- function(spots, pat_type, high_exp_prop = 0.3, 
                                grad_percent = 0.3, gene_id = NULL, seed =NULL){ 
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  x1 <- min(spots$x) ; x2 <- max(spots$x)
  y1 <- min(spots$y) ; y2 <- max(spots$y)
  
  # Calculate dynamic ranges
  x_range <- x2 - x1 ; y_range <- y2 - y1
  
  # Normalize coordinates to [0,1] for scale-invariant patterns
  spots$x_norm <- (spots$x - x1) / x_range
  spots$y_norm <- (spots$y - y1) / y_range
  
  # Determine number of high-expression spots for the pattern
  num_high <- round(high_exp_prop * nrow(spots))
  
  spots$eta <- 0  # Default to low expression
  
  if (is.null(pat_type)) {
    spots$eta <- 0
  } else {
    
    # Center in normalized space
    center_x <- 0.5  # Center in [0,1] normalized space
    center_y <- 0.5
    
    # Gene-specific random centers in normalized space
    if (!is.null(gene_id)) {
      center_x <- runif(1, 0.1, 0.9)  # Random center in [0.1, 0.9]
      center_y <- runif(1, 0.1, 0.9)
    }
    
    if (pat_type == "hotspot") {
      # Use normalized coordinates for distance calculation
      radius <- sqrt(sort((spots$x_norm - center_x)^2 + (spots$y_norm - center_y)^2)[num_high])
      spots$eta <- as.integer((spots$x_norm - center_x)^2 + (spots$y_norm - center_y)^2 <= radius^2)
      
    } else if (pat_type == "streak") {
      streak_x <- if (is.null(gene_id)) center_x else runif(1, 0.1, 0.9)
      streak_spots <- order(abs(spots$x_norm - streak_x))[1:num_high]  # Select num_high closest spots
      spots$eta[streak_spots] <- 1
      
    } else if (pat_type == "gradient") {
      threshold_x <- grad_percent  # In normalized space [0,1]
      spots$eta <- pmax(0, pmin(1, (spots$x_norm - threshold_x) / (1 - threshold_x)))
      if (!is.null(gene_id)) {
        shift <- runif(1, -0.2, 0.2)  # Smaller shift to avoid collapsing
        spots$eta <- pmax(0, pmin(1, spots$eta + shift))
      }
      
    } else if (pat_type == "ring") {
      distances <- (spots$x_norm - center_x)^2 + (spots$y_norm - center_y)^2
      ring_spots <- order(abs(distances - median(distances)))[1:num_high]
      spots$eta[ring_spots] <- 1
      
    } else if (pat_type == "wave") {
      phase <- if (!is.null(gene_id)) runif(1, -1, 1) else 0
      freq <- if (!is.null(gene_id)) runif(1, 2.5, 5.5) else 4
      w_y <- sin(freq * pi * (spots$x_norm + phase))
      # Scale wave to normalized y-space
      y_scaled <- 0 + (w_y - min(w_y)) / (max(w_y) - min(w_y)) * 1
      dist_to_curve <- abs(spots$y_norm - y_scaled)
      wave_spots <- order(dist_to_curve)[1:num_high]  # Select num_high closest spots
      spots$eta[wave_spots] <- 1
    }
  }
  # Assign expression regions based on pattern type
  if (is.null(pat_type)) {
    spots$mark <- "low"
  } else if (pat_type == "gradient") {
    spots$mark <- ifelse(spots$eta >= 0.75, "high",
                         ifelse(spots$eta <= 0.40, "low", "medium"))
  } else {
    spots$mark <- ifelse(spots$eta >= 0.75, "high", "low")
  }
  # Remove normalized columns to keep output clean
  spots$x_norm <- NULL ;  spots$y_norm <- NULL
  
  return(spots) # Return only pattern data
}



# ------------------------------------------------------------
# Function -- Pattern-based ST count data generator
# ------------------------------------------------------------

generate_genedata_pattern <- function( 
    spots, num_genes, prop = c(0.8, 0.2), pat_type, seed = NULL, G, gene_grp, cor_grp, tau_beta,
    high_exp_prop = 0.3,       # proportion of spots in high-expression region
    grad_percent  = 0.2,       # gradient threshold (for gradient pattern)
    beta0         = 2,         # log-scale baseline (Non-SE and low SE spots)
    beta1         = NULL,      # log-scale high-expression for SE genes
    sigma         = 0.3) {     # noise sd on log scale
  
  if (!is.null(seed)) set.seed(seed)
  
  if (length(prop) != 2 || any(prop < 0) || abs(sum(prop) - 1) > 1e-6)
    stop("prop must be a vector of length 2 summing to 1")
  
  spots <- as.data.frame(spots)
  colnames(spots)[1:2] <- c("x", "y")
  
  n <- nrow(spots)
  p <- num_genes
  
  ## default beta1: 3x fold change over beta0
  if (is.null(beta1)) beta1 <- beta0 + log(3)
  
  ## -----------------------------
  ## 1. Gene counts per category
  ## -----------------------------
  num_cat <- floor(p * prop)
  num_cat[which.max(prop * p - num_cat)] <- num_cat[which.max(prop * p - num_cat)] +
    (p - sum(num_cat))
  n_nonse <- num_cat[1]
  n_se    <- num_cat[2]
  
  ## -----------------------------
  ## 2. Random gene arrangement
  ## -----------------------------
  se_idx    <- sort(sample(seq_len(p), n_se))
  nonse_idx <- setdiff(seq_len(p), se_idx)
  Z         <- rep(1L, p)
  Z[se_idx] <- 2L
  
  ## -----------------------------
  ## 3. Gene-level dependence (CAR)
  ##    Scale Beta to [-0.3, 0.3] so it
  ##    acts as a small gene-level offset
  ##    without exploding log-lambda
  ## -----------------------------
  Beta <- create_car_beta(p, gene_grp, cor_grp, tau_beta)
  
  ## -----------------------------
  ## 4. Noise epsilon
  ## -----------------------------
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)
  
  ## -----------------------------
  ## 5. Generate counts
  ## -----------------------------
  result_mat <- matrix(0L, nrow = n, ncol = p)
  
  # Non-SE genes: flat Poisson with beta0
  for (j in nonse_idx) {
    loglambda_j     <- beta0 + Beta[j] + epsilon[, j]
    result_mat[, j] <- rpois(n, lambda = exp(loglambda_j))
  }
  
  # SE genes: Poisson with beta1 in high spots, beta0 elsewhere
  for (j in se_idx) {
    pattern_data <- simulate_expression(
      spots,  pat_type = pat_type,high_exp_prop = high_exp_prop,
      grad_percent  = grad_percent, gene_id  = j)
    
    loglambda_j              <- beta0         + Beta[j] + epsilon[, j]
    high_mask                <- pattern_data$mark == "high"
    med_mask                 <- pattern_data$mark == "medium"
    loglambda_j[high_mask]   <- beta1             + Beta[j] + epsilon[high_mask, j]
    loglambda_j[med_mask]    <- (beta0 + beta1)/2 + Beta[j] + epsilon[med_mask,  j]
    
    result_mat[, j] <- rpois(n, lambda = exp(loglambda_j))
  }
  
  ## -----------------------------
  ## 6. Gene naming NonSE: Gene{i}_NonSE SE:    Gene{i}_SE_{pat_type}
  ## -----------------------------
  gene_names            <- character(p)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx]    <- paste0("Gene", se_idx,    "_SE_", pat_type)
  colnames(result_mat)  <- gene_names
  
  ## -----------------------------
  ## Summary
  ## -----------------------------
  cat("---- Pattern Simulation Summary ----\n")
  cat("n spots   :", n, " | P genes:", p, "\n")
  cat("SE genes  :", n_se, " | Non-SE genes:", n_nonse, "\n")
  cat("SE indices:", paste(se_idx, collapse = ", "), "\n")
  cat("Pattern   :", pat_type, "\n")
  cat("beta0 =", round(beta0, 3), "| beta1 =", round(beta1, 3),
      "| fold change =", round(exp(beta1 - beta0), 2), "x\n")
  cat("sigma     :", sigma, "\n")
  cat("Observed counts: median =", median(result_mat),
      "| 95% =", paste(round(quantile(as.vector(result_mat), c(.025, .975))), collapse = ", "),
      "| max =", max(result_mat), "\n")
  cat("------------------------------------\n")
  
  return(list(
    Y         = result_mat,   # N x P count matrix
    se_idx    = se_idx,       # column indices of SE genes
    nonse_idx = nonse_idx,    # column indices of Non-SE genes
    Beta      = Beta,         # length P CAR pathway effects
    epsilon   = epsilon,      # N x P noise matrix
    beta0     = beta0,        # baseline log-mean
    beta1     = beta1         # SE high-region log-mean
  ))
}



# ------------------------------------------------------------
# Function -- SPHERE model data generator 
# ------------------------------------------------------------
gen_genedata_model <- function(
    spots, num_genes, prop = c(0.8, 0.2), G, gene_grp, 
    rho = 0.85, tau_beta = 10, tau_gs = NULL, ell_gs = NULL, 
    nugget = 1e-6, eps_sd = 0.20, 
    depth_model = c("poisson", "negbin", "fixed"), 
    depth_lambda = 5,    # mean for poisson or negbin
    depth_size   = 5,    # overdispersion for negbin (smaller = more overdispersed)
    depth_fixed  = 1, 
    target_mean  = 10,   # desired median count -- only parameter you may want to tune
    seed = NULL, verbose = TRUE) {
  
  depth_model <- match.arg(depth_model)
  
  # -----------------------
  # Checks
  # -----------------------
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
  
  # -----------------------
  # Sequencing depth N_i
  # -----------------------
  if (depth_model == "poisson") {
    N_i <- rpois(n, lambda = depth_lambda) + 1L
  } else if (depth_model == "negbin") {
    # Negative Binomial: mean = depth_lambda, variance = depth_lambda + depth_lambda^2/depth_size
    # depth_size controls overdispersion: smaller = more overdispersed
    N_i <- rnbinom(n, mu = depth_lambda, size = depth_size) + 1L
  } else {
    N_i <- rep.int(as.integer(depth_fixed), n)
  }
  
  # -------------------------------------------
  # Gene-level state Z_j : 1 = Non-SE, 2 = SE
  # -------------------------------------------
  n_se  <- round(P * prop[2])
  se_idx <- sort(sample(seq_len(P), n_se))
  nonse_idx <- setdiff(seq_len(p), se_idx)
  Z      <- rep(1L, P)
  Z[se_idx] <- 2L
  
  
  # -----------------------
  # CAR Beta (pathway effect)
  # -----------------------
  Beta <- create_car_beta2(num_genes = P, gene_grp = gene_grp, 
                           rho = rho, tau_beta = tau_beta)
  
  # -----------------------
  # iid noise epsilon_ij
  # -----------------------
  epsilon <- matrix(rnorm(n * P, mean = 0, sd = eps_sd), nrow = n, ncol = P)
  
  # ---------------------------------------------------
  # Spatial effect eta_ij only for SE genes, else 0
  # ---------------------------------------------------
  eta <- matrix(0, nrow = n, ncol = P)
  
  if (n_se > 0) {
    if (is.null(tau_gs)) tau_gs <- runif(n_se, 0.05, 0.30)
    if (is.null(ell_gs)) ell_gs <- runif(n_se, 0.50, 2.00)
    
    if (length(tau_gs) != n_se) stop("tau_gs must have length = number of SE genes.")
    if (length(ell_gs) != n_se) stop("ell_gs must have length = number of SE genes.")
    
    eta_se <- sim_gs_eta2(spots = spots[, 1:2, drop = FALSE], num_genes = n_se,
                          tau_gs = tau_gs, ell_gs = ell_gs, nugget = nugget)
    eta[, se_idx] <- eta_se
  } else {
    tau_gs <- numeric(0)
    ell_gs <- numeric(0)
  }
  
  # -----------------------
  # AUTO-CALIBRATE mu0
  # -----------------------
  # log lambda_ij = mu0 + Beta_j + epsilon_ij + eta_ij
  # We want: median(N_i * exp(mu0 + log_part)) = target_mean
  # => mu0 = log(target_mean) - log(median(N_i)) - median(log_part)
  # This is fully automatic -- no manual tuning needed when n or P changes
  
  log_part   <- sweep(epsilon + eta, 2, Beta, `+`)   # n x P
  mu0        <- log(target_mean) - log(median(N_i)) - median(log_part)
  
  # -----------------------
  # Generate counts
  # -----------------------
  lambda <- exp(mu0 + log_part)                        # latent rate n x P
  Y_mean <- outer(N_i, rep(1, P)) * lambda             # Poisson mean n x P
  Y      <- matrix(rpois(n * P, lambda = as.vector(Y_mean)), nrow = n, ncol = P)
  
  # Gene naming
  gene_names          <- character(P)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx]    <- paste0("Gene", se_idx,    "_SE_model")
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
        " | 95% =", paste(signif(quantile(as.vector(Y_mean), c(.025, .975)), 4), collapse = ", "),
        " | max =", signif(max(Y_mean), 4), "\n")
    cat("Observed counts      : median =", median(Y),
        " | 95% =", paste(quantile(as.vector(Y), c(.025, .975)), collapse = ", "),
        " | max =", max(Y), "\n")
    cat("-------------------------------\n")
  }
  
  return(list(
    Y       = Y,
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
## Plot all Figures
## ------------------------------------------------------------

pal <- colorRampPalette(c('#00274c', '#00274c',"lightyellow2",'#ffcb05','#ffcb05'))
# Function to create a ggplot for each element in the list and arrange them in a grid
plot_gene_data <- function(data, nrow = 1) {
  plots <- lapply(names(data), function(gene_name) {
    ggplot(data[[gene_name]], aes(x, y, color = expression)) +
      geom_point(size = 4) +
      scale_color_gradientn(colours = pal(5), limits = c(0, 1)) + 
      theme_minimal() +
      labs(title = paste("Relative expression Counts for \n", gene_name),
           color = "Relative \nexpression") +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  })
  
  # Arrange all plots in a grid with the specified number of rows
  do.call(grid.arrange, c(plots, nrow = nrow))
}



##---------------------------------------------------------
## Function to filter genes with less total counts 
##---------------------------------------------------------

filter_columns_by_sum <- function(df, threshold = 10, keep_first = TRUE) {
  # Convert to data.table for convenience
  library(data.table)
  dt <- as.data.table(df)
  
  # Identify columns to check (skip first if keep_first = TRUE)
  if (keep_first) {
    numeric_part <- dt[, -1, with = FALSE]
    sums <- colSums(numeric_part, na.rm = TRUE)
    cols_keep <- c(names(dt)[1], names(numeric_part)[sums > threshold])
  } else {
    sums <- colSums(dt, na.rm = TRUE)
    cols_keep <- names(dt)[sums > threshold]
  }
  # Return filtered table
  return(dt[, ..cols_keep])
}

##---------------------------------------------------------
## Function to make RBF basis matrix
##---------------------------------------------------------
make_rbf_basis <- function(coords, r = 30, lengthscale = 1.5) {
  N <- nrow(coords) # choose r knot centers using k-means or random subset
  km <- kmeans(coords, centers = r)
  knots <- km$centers
  # compute pairwise distance matrix between points and knots
  d <- as.matrix(dist(rbind(coords, knots)))
  D <- d[1:N, (N+1):(N+r)]  # N x r
  # set lengthscale = median distance between knots if not provided
  if (is.null(lengthscale)) {
    lengthscale <- median(as.matrix(dist(knots)))  }
  # RBF basis
  Phi <- exp( - (D^2) / (2 * lengthscale^2) )
  return(Phi)   # N x r
}

##---------------------------------------------------------
## Function to make RBF distance matrix
##---------------------------------------------------------
make_rbf_dist <- function(coords, r = 30) {
  coords <- as.matrix(coords)
  N <- nrow(coords)
  km <- kmeans(coords, centers = r)
  knots <- km$centers
  all_d <- as.matrix(dist(rbind(coords, knots)))
  D <- all_d[1:N, (N+1):(N+r)]
  D
}


##---------------------------------------------------------
## Function to generate gene pathway group
##---------------------------------------------------------
generate_gene_grp <- function(p, G = 3) {
  reps <- ceiling(p / G)
  rep(1:G, each = reps)[1:p]
}


##---------------------------------------------------------
## Function to get the tau and ell for SE genes
##---------------------------------------------------------
get_tau_ell <- function(p, prop) {
  se_count <- round(prop[2] * p)
  list(
    tau_gs = t_gs[1:se_count],
    ell_gs = l_gs[1:se_count]
  )
}

##---------------------------------------------------------
## relative expression
##---------------------------------------------------------

relative_expr <- function(raw_exp) {
  raw_exp_log <- log1p(raw_exp) # log(1 + x)
  (raw_exp_log - min(raw_exp_log)) / (max(raw_exp_log) - min(raw_exp_log))
}
