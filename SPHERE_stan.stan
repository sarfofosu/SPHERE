// ==================================================================================
// SPHERE: A Spatial Poisson Hierarchical modEl with pathway-infoRmed gEne networks
// ==================================================================================
// This Stan model identifies spatially expressed genes (SVGs) in spatial
// transcriptomics data. It uses:
//   1. A Poisson log-normal observation model for count data
//   2. A low-rank Gaussian Process (GP) to capture spatial patterns
//   3. A gene-level two-component mixture to classify genes as
//      spatially expressed (SE) or non-spatially expressed (non-SE)
//   4. A Conditional Autoregressive (CAR) prior that shares information
//      across genes belonging to the same biological pathway
// =============================================================================



// ---------------------------------------------------------------
// DATA BLOCK
// ---------------------------------------------------------------
// Everything here is fixed and supplied by the user from R.
// Stan reads these values once before sampling begins.
// ---------------------------------------------------------------
data {
  // === Core Dimensions ===
  int<lower=1> N;                              // number of spots
  int<lower=1> P;                              // number of genes
  array[N,P] int<lower=0> Y;                   // Gene expression count matrix (spots x genes)
  vector<lower=0>[N] N_i;                      // Normalization offsets : Accounts for differences in sequencing depth across spots

  // === Low-Rank GP Basis Function Inputs ===
  // Instead of computing a full N x N covariance matrix (expensive),
  // we approximate the GP using r << N inducing knots.
  
  int<lower=1> r;                              // Number of knots (inducing points) for the low-rank GP, r << N
  matrix[N,r] D;                               // Squared Euclidean distances from each spot to each knot
  matrix[N,r] Phi;                             // Basis function matrix (N spots x r knots)  Phi[i,b] = basis function value linking spot i to knot b

  // === Hyperparameters for Prior Distributions ===
  // These are fixed constants that control the shape of the priors.
  
  real<lower=0> a_err;  real<lower=0> b_err;    // Location (mean) and Scale (sd) for the half-normal prior on sigma_sd[j]
  vector<lower=0>[2] alpha;                     // Dirichlet concentration parameters (length 2) alpha = (alpha_1, alpha_2) and Controls the prior on mixture weights pii[j]
  real<lower=0> a_gs;   real<lower=0> b_gs;     // Location (mean) and Scale (sd) for the half-normal prior on sig_eta_gs[j]
  real<lower=0> a_gsl;  real<lower=0> b_gsl;    // Mean and sd parameter for log-normal prior on ell_gs[j] (GP lengthscale)
  real<lower=0> a_mu;  real<lower=0> b_mu;      // Mean and sd for the normal prior on mu0 (global intercept)

  // === Pathway / CAR Grouping ===
  // Genes are organized into G biological pathway groups.
  // Genes in the same group are treated as "neighbors" in the CAR prior,
  // allowing them to share information about expression effects.
  int<lower=1> G;                                            // Number of distinct gene groups (pathways) 
  array[P] int<lower=1, upper=G> gene_group;                 // Pathway membership for each gene; gene_group[j] = which group gene j belongs to
  real<lower=0> a_tau_beta; real<lower=0> b_tau_beta;        // Location and scale for half-normal prior on sigma_beta (CAR precision)
  real<lower=0> a_rho;      real<lower=0> b_rho;             // shapes parameter for Beta prior; rho controls strength of genes in the CAR
}


// ---------------------------------------------------------------
// PARAMETER
// ---------------------------------------------------------------

parameters {
  real mu0;
  // gene-level mixture probabilities
  array[P] simplex[2] pii;
  // latent log-rates (this is log(lambda_ij) in your theory)
  matrix[N, P] loglambda;
  // observation-level noise (per gene)
  vector<lower=1e-6>[P] sigma_sd;
  // spatial GP amplitude + lengthscale (per gene)
  vector<lower=1e-6>[P] sig_eta_gs;
  vector<lower=1e-6>[P] ell_gs;
  // low-rank GP weights (per gene)
  matrix[r, P] w;
  // pathway CAR gene effects
  vector[P] Beta;
  real<lower=1e-6> sigma_beta;
  real<lower=1e-6, upper=0.999> rho;
}



transformed parameters {
  //-----------------------------------------------------------
  // Compute low-rank RBF spatial effects: use median of distance as lengthscale
  //   eta_gs_trans[:, j] = sig_eta_gs[j] * (Phi * eta_gs[:, j])
  //-----------------------------------------------------------
    matrix[N, P] eta;
    matrix[N, P] tmp = Phi * w;  // N x P
    for (j in 1:P) {
      eta[, j] = sig_eta_gs[j] * tmp[, j];
    }
}



model {
  // --------------------------
  // Priors
  // --------------------------
  target += normal_lpdf(mu0 | a_mu, b_mu);

  for (j in 1:P) {
    pii[j] ~ dirichlet(alpha);

    // “half-normal” via normal with lower bound
    sigma_sd[j] ~ normal(a_err, b_err);
    sig_eta_gs[j] ~ normal(a_gs, b_gs);
    ell_gs[j] ~ lognormal(a_gsl, b_gsl);

    // GP weights
    w[, j] ~ normal(0, 1);
  }

  sigma_beta ~ normal(a_tau_beta, b_tau_beta);
  rho ~ beta(a_rho, b_rho);

  // --------------------------
  // Conditional CAR prior for Beta (pathway neighbors)
  // Matches: Beta_j | neighbors ~ N( rho * mean(neighbors), sigma_beta / sqrt(n_j) )
  // --------------------------
  {
    // adjacency is implicit: same gene_group => neighbors
    for (j in 1:P) {
      real nj = 0;
      real sum_nb = 0;

      for (k in 1:P) {
        if (k != j && gene_group[k] == gene_group[j]) {
          nj += 1;
          sum_nb += Beta[k];
        }
      }

      if (nj > 0.5) {
        Beta[j] ~ normal(rho * (sum_nb / nj), sigma_beta / sqrt(nj));
      } else {
        Beta[j] ~ normal(0, 10); // isolated gene
      }
    }
  }

  // --------------------------
  // Mixture PRIOR for loglambda (GENE-LEVEL, not spot-level)
  // This is the critical fix: Z_j is per gene.
  //
  // Non-SE component (Z_j=1): loglambda_ij ~ Normal(mu0 + Beta_j, sigma_sd_j)
  // SE component (Z_j=2):     loglambda_ij ~ Normal(mu0 + Beta_j + eta_ij, sigma_sd_j)
  //
  // We marginalize Z_j via log_mix over the WHOLE gene likelihood (sum over i).
  // --------------------------

    for (j in 1:P) {
    for (i in 1:N) {
      target += poisson_log_lpmf(Y[i,j] |log(N_i[i]) + loglambda[i,j]);
      target += log_mix(
                    pii[j,1], normal_lpdf(loglambda[i,j] | mu0 + Beta[j],             sigma_sd[j]),
                              normal_lpdf(loglambda[i,j] | mu0 + Beta[j] + eta[i, j], sigma_sd[j]));
    }
  }
}


// ---------------------------------------------------------------
// GENERATED QUANTITIES BLOCK - POSTERIOR CLASSIFICATION
// ---------------------------------------------------------------
// Computed AFTER each HMC draw, using the sampled parameter values.
// These are derived quantities for posterior inference 
// ---------------------------------------------------------------
generated quantities {
  array[P] int<lower=1, upper=2> Z;                      // Posterior gene-level classification indicator: Z[j] = 1 NON-SE ; Z[j] = 2 -> SE
  for (j in 1:P) {
    vector[2] lpp = rep_vector(0,2);                     // Initialize log-posterior-predictive scores for each component
    for (i in 1:N) {
      lpp[1] += normal_lpdf(loglambda[i,j] | mu0 + Beta[j], sigma_sd[j]);
      lpp[2] += normal_lpdf(loglambda[i,j] | mu0 + Beta[j] +  eta[i,j], sigma_sd[j]);
    }
    lpp += log(pii[j]);
    Z[j] = categorical_logit_rng(lpp);      // Over many HMC samples, the fraction of times Z[j] = 2
                                            // gives the posterior probability that gene j is an SVG.
  }
}

