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

  real<lower=0> mu_noise; real<lower=0> sd_noise;                   // scale of half-normal prior on sigma_sd
  vector<lower=0>[2] alpha;                                         // Dirichlet concentration parameters (length 2) alpha = (alpha_1, alpha_2) and Controls the prior on mixture weights pii[j]
  real<lower=0> mu_gp_amplitude; real<lower=0> sd_gp_amplitude;     // mean and scale of half-normal prior on sig_eta_gs[j]
  real<lower=0> mu_gp_lengthscale; real<lower=0> sd_gp_lengthscale; // mean and scale of log-normal prior on ell_gs
  real mu_intercept; real<lower=0> sd_intercept;                    // scale of normal prior on mu0

  // === Pathway / CAR Grouping ===
  // Genes are organized into G biological pathway groups.
  // Genes in the same group are treated as "neighbors" in the CAR prior,
  // allowing them to share information about expression effects.
  int<lower=1> G;                                            // Number of distinct gene groups (pathways)
  array[P] int<lower=1, upper=G> gene_group;                 // Pathway membership for each gene; gene_group[j] = which group gene j belongs to
  real<lower=0> mu_car_precision; real<lower=0> sd_car_precision; // mean and scale of half-normal prior on sigma_beta
  real<lower=0> shape_car_cor; real<lower=0> rate_car_cor;        // first and second shape of Beta prior on rho_CAR
}


// ---------------------------------------------------------------
// PARAMETERS BLOCK
// ---------------------------------------------------------------
// These are the unknown quantities that Stan will estimate via HMC.
// Stan explores the joint posterior distribution over all of these.
// ---------------------------------------------------------------

parameters {
  real mu0;                               // Global mean log-expression intercept; Shared across all genes and all spots; Represents the baseline expression level on the log scale
  array[P] simplex[2] pii;                // Per-gene mixture probabilities pii[j,1] = prob gene j is NON-SE (Z_j=1), pii[j,2] = prob gene j IS SE (Z_j=2)
  matrix[N, P] loglambda;                 // Latent log-rate matrix (spots x genes), where lambda_ij is the true expression rate for gene j at spot i
  vector<lower=1e-6>[P] sigma_sd;         // Per-gene observation-level standard deviation; Controls how much loglambda[i,j] can deviate from its mean
  vector<lower=1e-6>[P] sig_eta_gs;       // Per-gene GP amplitude (signal standard deviation)  Scales the magnitude of the spatial effect for gene j
  vector<lower=1e-6>[P] ell_gs;           // Per-gene GP lengthscale; Controls how quickly spatial correlation decays with distance
  matrix[r, P] w;                         // Low-rank GP weight matrix (knots x genes); w[b,j] = weight for knot b contributing to gene j's spatial effect
  vector[P] Beta;                         // Gene-level pathway effect coefficients; Beta[j] = the effect of gene j's pathway membership on expression
  real<lower=1e-6> sigma_beta;            // CAR conditional standard deviation; Controls how tightly genes within a pathway are pulled together
  real<lower=1e-6, upper=0.999> rho;      // Spatial autocorrelation parameter for the CAR prior; Controls the degree of borrowing between pathway neighbors
}



// ---------------------------------------------------------------
// TRANSFORMED PARAMETERS BLOCK
// ---------------------------------------------------------------
// Deterministic transformations of parameters, computed each iteration.
// These are intermediate quantities derived from the raw parameters.
// ---------------------------------------------------------------

transformed parameters {
  matrix[N, P] eta;                            // The spatial effect matrix (spots x genes); eta[i,j] = spatial contribution to log-rate for gene j at spot i
  matrix[N, P] tmp = Phi * w;                  // Intermediate: raw basis expansion (N x P); tmp = Phi (N x r) times w (r x P) = N x P matrix; use median of distance as lengthscale
  for (j in 1:P) {
    eta[, j] = sig_eta_gs[j] * tmp[, j];       // Scale the raw spatial effect by the gene-specific GP amplitude; This improves HMC sampling by decoupling the amplitude
  }
}




// ---------------------------------------------------------------
// MODEL BLOCK
// ---------------------------------------------------------------
// Specifies the full joint log-posterior:
// log p(parameters | data) = log p(data | parameters) + log p(parameters)
// Stan adds all target += statements to compute this.
// ---------------------------------------------------------------

model {
  // ===========================================
  // SECTION 1: PRIOR DISTRIBUTIONS
  // ===========================================
  // These encode our prior beliefs about parameter values before seeing data.

  target += normal_lpdf(mu0 | mu_intercept, sd_intercept);        // Prior on global intercept: mu0 ~ Normal(a_mu, b_mu)
  for (j in 1:P) {                                                // Loop over all P genes
    pii[j] ~ dirichlet(alpha);                                    // Prior on mixture weights for gene j; pii[j] ~ Dirichlet(alpha[1], alpha[2]) and Controls the prior probability of being spatially variable
    sigma_sd[j] ~ normal(mu_noise, sd_noise);                     // Prior on observation noise for gene j; Larger b_err = more diffuse prior, allows more noise
    sig_eta_gs[j] ~ normal(mu_gp_amplitude, sd_gp_amplitude);     // Prior on GP amplitude for gene j; Controls expected strength of spatial signal
    ell_gs[j] ~ lognormal(mu_gp_lengthscale, sd_gp_lengthscale);  // Prior on GP lengthscale for gene j
    w[, j] ~ normal(0, 1);                                        // Prior on GP weights for gene j;  Standard normal prior on each of the r knot weights
  }
  sigma_beta ~ normal(mu_car_precision, sd_car_precision);        // Prior on CAR conditional SD;  Controls how dispersed gene effects are within pathways
  rho ~ beta(shape_car_cor, rate_car_cor);                        // Prior on CAR autocorrelation parameter; Larger a_rho pushes rho toward 1 (stronger spatial smoothing)


  // ================================================
  // SECTION 2: CAR PRIOR FOR GENE PATHWAY EFFECTS
  // ================================================
  // This implements a Conditional Autoregressive (CAR) prior on Beta[j].
  // Genes in the same pathway are "neighbors" and share information.
  // The full conditional for each gene is: Beta[j] | Beta[-j] ~ Normal(rho * mean(neighbors), sigma_beta / sqrt(n_j)), where n_j = number of neighbors of gene j.
  {
    for (j in 1:P) {                                  // Loop over each gene
      real nj = 0;                                    // Counter: number of neighbors of gene j
      real sum_nb = 0;                                // Accumulator: sum of Beta values of gene j's neighbors

      for (k in 1:P) {                                // Inner loop: check all other genes
        if (k != j && gene_group[k] == gene_group[j]) {// Gene k is a neighbor of gene j if:
          nj += 1;                                    // Increment neighbor count
          sum_nb += Beta[k];                          // Add gene k's Beta to the running sum
        }
      }
      if (nj > 0.5) {                                 // If gene j has at least 1 neighbor (using > 0.5 instead)
        Beta[j] ~ normal(rho * (sum_nb / nj), sigma_beta / sqrt(nj));
                                                      // CAR full conditional
      } else {
        Beta[j] ~ normal(0, 10);                      // Isolated gene (no pathway neighbors) Given a weakly informative prior centered at 0
      }
    }
  }


  // =============================================================
  // SECTION 3: MIXTURE LIKELIHOOD
  // =============================================================
  // This is the core of the model. For each gene j and spot i:
  //   1. The observed count Y[i,j] follows a Poisson with rate N_i * lambda_ij
  //   2. The latent log-rate loglambda[i,j] comes from a two-component mixture:
  //      - Component 1 (non-SE): loglambda ~ Normal(mu0 + Beta[j], sigma_sd[j])
  //      - Component 2 (SE):     loglambda ~ Normal(mu0 + Beta[j] + eta[i,j], sigma_sd[j])
  for (j in 1:P) {               // Loop over all genes
    for (i in 1:N) {             // Loop over all spots

      target += poisson_log_lpmf(Y[i,j] | log(N_i[i]) + loglambda[i,j]); // Y[i,j] ~ PoissonLog(log(N_i[i]) + loglambda[i,j])
      target += log_mix(
                    pii[j,1],                                            // MIXTURE PRIOR on loglambda (marginalized over Z_j):
                    normal_lpdf(loglambda[i,j] | mu0 + Beta[j],             sigma_sd[j]),  //         -> mean is just baseline + pathway effect
                    normal_lpdf(loglambda[i,j] | mu0 + Beta[j] + eta[i, j], sigma_sd[j])); //         -> mean shifts by the GP-derived spatial pattern
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

