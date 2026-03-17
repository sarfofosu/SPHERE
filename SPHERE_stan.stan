// ================================================================
//  STAN model for SPHERE
// ================================================================



// ---------------------------------------------------------------
// DATA
// ---------------------------------------------------------------
data {
  // Data Setup
  int<lower=1> N;                              // number of spots
  int<lower=1> P;                              // number of genes
  array[N,P] int<lower=0> Y;                  // gene expression counts
  vector<lower=0>[N] N_i;                      // Normalization offsets

  // RBF Low-rank GP basis
  int<lower=1> r;                              // number of knots / basis functions
  matrix[N,r] D;                               // squared distances ||s_i - k_b||^2  (N x r)
  matrix[N,r] Phi;                             // distances from spots to knots

  // Hyperparameters for Priors
  real<lower=0> a_err;  real<lower=0> b_err;    // sigma_sd
  vector<lower=0>[2] alpha;                     // Dirichlet for pi_j
  real<lower=0> a_gs;   real<lower=0> b_gs;     // sig_eta_gs
  real<lower=0> a_gsl;  real<lower=0> b_gsl;    // ell_gs
  real<lower=0> a_mu;  real<lower=0> b_mu;      // mu0

  // Pathway / CAR grouping
  int<lower=1> G;
  array[P] int<lower=1, upper=G> gene_group;
  real<lower=0> a_tau_beta; real<lower=0> b_tau_beta;
  real<lower=0> a_rho;      real<lower=0> b_rho;
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
// POSTERIOR CLASSIFICATION
// ---------------------------------------------------------------
generated quantities {
  array[P] int<lower=1, upper=2> Z;
  for (j in 1:P) {
    vector[2] lpp = rep_vector(0,2);
    for (i in 1:N) {
      lpp[1] += normal_lpdf(loglambda[i,j] | mu0 + Beta[j], sigma_sd[j]);
      lpp[2] += normal_lpdf(loglambda[i,j] | mu0 + Beta[j] +  eta[i,j], sigma_sd[j]);
    }
    lpp += log(pii[j]);
    Z[j] = categorical_logit_rng(lpp);
  }
}

