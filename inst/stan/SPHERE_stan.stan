// ==================================================================================
// SPHERE: A Spatial Poisson Hierarchical modEl with pathway-infoRmed gEne networks
// ==================================================================================
// This Stan model identifies spatially expressed genes (SVGs) in spatial
// transcriptomics data. It uses:
//   1. A Poisson log-normal observation model for count data
//   2. A low-rank Gaussian Process (GP) to capture spatial patterns
//   3. A gene-level Bernoulli model to classify genes as
//      spatially expressed (SE) or non-spatially expressed (non-SE)
//   4. A Conditional Autoregressive (CAR) prior (Leroux or sparse proper) that
//      shares information on the SE log-odds across genes belonging to the
//      same biological pathway
//
// MODEL OPTIONS:
//   use_shared_Phi = 1: Phi is supplied directly from R and shared by all genes;
//                       the GP basis weights w are integrated out analytically
//   use_shared_Phi = 0: construct Phi_j using D and gene-specific ell_gs[j]
//                       (basis weights w sampled explicitly)
//   use_pathway  = 1: pathway-informed detection;  0: no pathway          (M3)
//   fix_rho_zero = 1: rho_u = 0  (iid u, pathway-free)                    (M2)
//   fix_rho_zero = 0: estimate rho_u                                      (M1)
//   use_car_lx   = 1: Leroux CAR on u_j
//   use_car_lx   = 0: sparse proper CAR on u_j
// ==================================================================================



// ---------------------------------------------------------------
// FUNCTIONS BLOCK
// ---------------------------------------------------------------
// User-defined log densities for the two pathway CAR priors.
// Both return the exact log density up to an additive constant and
// use a sparse edge list (W_sparse) so W is never stored as a full matrix.
// ---------------------------------------------------------------
functions {
  // Leroux CAR: u ~ N(0, Q^{-1}),  Q = tau_u [ (1-rho_u) I + rho_u (D - W) ]
  // Exact log density up to the additive constant -(P_car/2) log(2 pi).
  real leroux_car_lpdf(
      vector u_car,real sigma_u,real rho_u,array[,] int W_sparse,
      vector D_sparse,vector lap_eig,int P_car,int W_n) {
    real tau_u=inv_square(sigma_u);                              // Precision parameter: tau_u = 1 / sigma_u^2
    vector[P_car] Wu=rep_vector(0.0,P_car);                      // W * u; accumulates the sum of each gene's neighbour effects
    for(e in 1:W_n){                                             // Loop over each edge (pair of neighbouring genes)
      int j1=W_sparse[e,1];                                      // First gene in the edge
      int j2=W_sparse[e,2];                                      // Second gene in the edge
      Wu[j1]+=u_car[j2];                                         // Add neighbour j2's effect to gene j1
      Wu[j2]+=u_car[j1];                                         // Add neighbour j1's effect to gene j2 (W is symmetric)
    }
    real log_det=0.0;                                            // log|(1-rho)I + rho(D-W)|
    for(j in 1:P_car){                                           // Log-determinant computed from Laplacian eigenvalues
      log_det+=log1p(rho_u*(lap_eig[j]-1));                      // log(1 - rho_u + rho_u * lap_eig[j]); avoids an O(P^3) determinant each iteration
    }
    return 0.5*(P_car*log(tau_u)+log_det                         // Normalising terms: log|Q| = P log(tau_u) + log_det
                -tau_u*((1-rho_u)*dot_self(u_car)                // Quadratic form u' Q u: independent part (1-rho_u) u'u
                        +rho_u*(dot_product(u_car.*D_sparse,u_car)-dot_product(u_car,Wu))));  // plus the graph part rho_u (u'Du - u'Wu)
  }

  // Sparse proper CAR: u ~ N(0, Q^{-1}),  Q = tau_u (D - rho_u W)
  real sparse_car_lpdf(
      vector u_car,real sigma_u,real rho_u,array[,] int W_sparse,
      vector D_sparse,vector car_eig,int P_car,int W_n) {
    real tau_u=inv_square(sigma_u);                              // Precision parameter
    vector[P_car] Wu=rep_vector(0.0,P_car);                      // W * u
    if(W_n>0){                                                   // Skip when the graph has no edges
      for(e in 1:W_n){                                           // Sparse multiplication W * u
        int j1=W_sparse[e,1];                                    // First gene in the edge
        int j2=W_sparse[e,2];                                    // Second gene in the edge
        Wu[j1]+=u_car[j2];                                       // Add neighbour j2's effect to gene j1
        Wu[j2]+=u_car[j1];                                       // Add neighbour j1's effect to gene j2 (W is symmetric)
      }
    }
    real log_det=0.0;                                            // log|D-rho_uW| contribution
    for(j in 1:P_car){                                           // Log-determinant computed from the eigenvalues of D^{-1/2} W D^{-1/2}
      log_det+=log1m(rho_u*car_eig[j]);                          // log(1 - rho_u * car_eig[j]); the constant sum(log D) is dropped
    }
    real u_D_u=dot_product(u_car.*D_sparse,u_car);               // u' D u
    real u_W_u=dot_product(u_car,Wu);                            // u' W u
    // Exact proper CAR log density, up to an additive constant
    return 0.5*(P_car*log(tau_u)+log_det-tau_u*(u_D_u-rho_u*u_W_u));
  }
}



// ---------------------------------------------------------------
// DATA BLOCK
// ---------------------------------------------------------------
// Everything here is fixed and supplied by the user from R.
// Stan reads these values once before sampling begins.
// ---------------------------------------------------------------
data {
  // === Core Dimensions ===
  int<lower=1> N;                                                // Number of spatial spots
  int<lower=1> P;                                                // Number of genes
  array[N,P] int<lower=0> Y;                                     // Gene expression count matrix (spots x genes)
  vector<lower=0>[N] N_i;                                        // Spot library sizes N(s_i) (> 0); accounts for differences in sequencing depth across spots

  // === Low-Rank GP Basis Function Inputs ===
  // Instead of computing a full N x N covariance matrix (expensive),
  // we approximate the GP using r << N inducing knots.
  int<lower=1> r;                                                // Number of knots (inducing points) for the low-rank GP, r << N
  matrix[N,r] Phi;                                               // Shared basis (uncentred) from R; Phi[i,b] links spot i to knot b
  matrix<lower=0>[N,r] D;                                        // Spot-to-knot distances (same knots); used to build gene-specific bases

  // === Model Options ===
  // Integer switches that select which version of the model is fitted.
  int<lower=0,upper=1> use_shared_Phi;                           // 1 = shared basis with w integrated out; 0 = gene-specific basis with w sampled
  int<lower=0,upper=1> use_pathway;                              // 1 = pathway-informed detection; 0 = no pathway (M3)
  int<lower=0,upper=1> fix_rho_zero;                             // 1 = rho_u fixed at 0 (M2); 0 = rho_u estimated (M1)
  int<lower=0,upper=1> use_car_lx;                               // 1 = Leroux CAR; 0 = sparse proper CAR
  int<lower=0> K_u;                                              // P when pathway used; 0 otherwise (sizes u)
  int<lower=0> K_sigma_u;                                        // 1 when pathway used; 0 otherwise (sizes sigma_u)
  int<lower=0> K_rho;                                            // 1 when rho_u estimated; 0 otherwise (sizes rho_u_free)

  // === Hyperparameters for Prior Distributions ===
  // These are fixed constants that control the shape of the priors.
  real<lower=0> mu_noise;                                        // Location of half-normal prior on sigma_sd
  real<lower=0> sd_noise;                                        // Scale of half-normal prior on sigma_sd
  real<lower=0> mu_gp_amplitude;                                 // Location of half-normal prior on sig_eta_gs
  real<lower=0> sd_gp_amplitude;                                 // Scale of half-normal prior on sig_eta_gs
  real mu_gp_lengthscale;                                        // Log-scale mean of log-normal prior on ell_gs
  real<lower=0> sd_gp_lengthscale;                               // Log-scale SD of log-normal prior on ell_gs
  real mu_intercept;                                             // Mean of normal prior on mu0
  real<lower=0> sd_intercept;                                    // SD of normal prior on mu0

  // === Pathway / CAR Structure ===
  // Genes in the same pathway are "neighbours" in the CAR graph.
  // Genes with no neighbours are handled separately as isolated genes.
  int<lower=0> P_car;                                            // Number of genes with at least one pathway neighbour
  array[P_car] int<lower=1,upper=P> car_index;                   // Indices of those genes among 1..P
  int<lower=0> P_iso;                                            // Number of isolated genes (no pathway neighbours)
  array[P_iso] int<lower=1,upper=P> iso_index;                   // Indices of isolated genes among 1..P
  int<lower=0> W_n;                                              // Number of edges in the pathway graph
  array[W_n,2] int<lower=1,upper=P_car> W_sparse;                // Edge list: each row is a pair of neighbouring genes
  vector<lower=1>[P_car] D_sparse;                               // Number of neighbours n_j for each connected gene
  vector[P_car] lap_eig;                                         // Laplacian eigenvalues of D - W (Leroux CAR)
  vector[P_car] car_eig;                                         // Eigenvalues of D^{-1/2} W D^{-1/2} (sparse CAR)

  // === Independent Beta Prior Hyperparameters ===
  real<lower=0> mu_beta_sig;                                     // Location of half-normal prior on sigma_beta
  real<lower=0> sd_beta_sig;                                     // Scale of half-normal prior on sigma_beta

  // === Pathway-Informed Detection Hyperparameters ===
  real alpha0_mean;                                              // Mean of normal prior on alpha0 (global SE log-odds)
  real<lower=0> alpha0_sd;                                       // SD of normal prior on alpha0
  real<lower=0> mu_u_sig;                                        // Location of half-normal prior on sigma_u
  real<lower=0> sd_u_sig;                                        // Scale of half-normal prior on sigma_u
  real<lower=0> shape_rho_u;                                     // rho_u ~ Beta(shape_rho_u, rate_rho_u)
  real<lower=0> rate_rho_u;                                      // Second shape of Beta prior on rho_u
}



// ---------------------------------------------------------------
// TRANSFORMED DATA BLOCK
// ---------------------------------------------------------------
// Quantities computed once from the data before sampling begins.
// These never change during MCMC, so computing them here saves time.
// ---------------------------------------------------------------
transformed data {
  vector[N] log_N_i=log(N_i);                                    // Log library sizes; used as the Poisson offset
  matrix[N,r] Phi_c;                                             // column-centred basis
  for(b in 1:r){                                                 // Centre each basis column so the GP carries only the spatial pattern
    Phi_c[,b]=Phi[,b]-mean(Phi[,b]);                             // and Beta[j] carries the overall expression level
  }
  // Spectral form of the shared basis for the collapsed SE likelihood.
  matrix[r,r] PtP=crossprod(Phi_c);                              // Phi_c' Phi_c (r x r)
  vector[r] lam_all=eigenvalues_sym(PtP);                        // All eigenvalues of Phi_c' Phi_c
  matrix[r,r] V_all=eigenvectors_sym(PtP);                       // Corresponding eigenvectors
  real lam_tol=1e-8*max(lam_all);                                // Tolerance for treating an eigenvalue as zero
  int r_eff=0;                                                   // Effective rank of the centred basis
  for(k in 1:r){
    if(lam_all[k]>lam_tol) r_eff+=1;                             // Count the non-negligible eigenvalues
  }
  vector[r_eff] lam;                                             // nonzero eigenvalues of Phi_c' Phi_c
  matrix[N,r_eff] G;                                             // orthonormal: Phi_c V Lam^{-1/2}
  {
    int k_eff=0;                                                 // Position in the reduced eigen-system
    for(k in 1:r){
      if(lam_all[k]>lam_tol){                                    // Keep only non-negligible directions
        k_eff+=1;
        lam[k_eff]=lam_all[k];                                   // Store the eigenvalue
        G[,k_eff]=Phi_c*V_all[,k]/sqrt(lam_all[k]);              // Store the matching orthonormal basis column
      }
    }
  }
  int K_w=(use_shared_Phi==1)?0:r;                               // explicit w only if gene-specific
}



// ---------------------------------------------------------------
// PARAMETERS BLOCK
// ---------------------------------------------------------------
// These are the unknown quantities that Stan will estimate via HMC.
// Stan explores the joint posterior distribution over all of these.
// ---------------------------------------------------------------
parameters {
  real mu0;                                                      // Global intercept; baseline log-expression shared across all genes and spots
  matrix[N,P] z_ll;                                              // standardized latent residuals; loglambda = mu0 + Beta_j + sigma_sd_j * z_ll
  vector<lower=1e-6>[P] sigma_sd;                                // Gene-specific residual SD; controls how far log-rates deviate from their mean
  vector<lower=1e-6>[P] sig_eta_gs;                              // Gene-specific GP amplitude; scales the strength of gene j's spatial effect
  vector<lower=1e-6>[P] ell_gs;                                  // Gene-specific GP lengthscales; controls how fast spatial correlation decays
  matrix[K_w,P] w;                                               // GP basis weights (gene-specific Phi only); empty when use_shared_Phi = 1
  vector[P] Beta;                                                // Independent gene-level effects; gene j's overall expression level
  real<lower=1e-6> sigma_beta;                                   // Gene-level effect SD; controls how spread out Beta is across genes
  real alpha0;                                                   // Global SE log-odds; sets the overall prevalence of SE genes
  vector[K_u] u;                                                 // Pathway-informed latent detection effects; shift each gene's SE log-odds
  vector<lower=1e-6>[K_sigma_u] sigma_u;                         // Empty when pathway is not used; SD of the pathway effects u
  vector<lower=1e-6,upper=0.999>[K_rho] rho_u_free;              // Empty unless rho_u estimated; strength of borrowing between pathway neighbours
}



// ---------------------------------------------------------------
// TRANSFORMED PARAMETERS BLOCK
// ---------------------------------------------------------------
// Deterministic transformations of parameters, computed each iteration.
// These are intermediate quantities derived from the raw parameters.
// ---------------------------------------------------------------
transformed parameters {
  vector<lower=0,upper=1>[P] pi;                                 // Gene-specific prior P(Z_j=1)
  vector[P] u_centered;                                          // Centered pathway effects
  real<lower=0,upper=0.999> rho_u;                               // CAR dependence for detection
  vector[P] logLR;                                               // log p(z_j|SE) - log p(z_j|non-SE)

  // === CAR dependence parameter ===
  if(use_pathway==0){
    rho_u=0;                                                     // No-pathway model (M3)
  }else if(fix_rho_zero==1){
    rho_u=0;                                                     // iid u (M2)
  }else{
    rho_u=rho_u_free[1];                                         // estimated (M1)
  }

  // === Gene-specific prior probability of being SE ===
  if(use_pathway==1){
    u_centered=u-mean(u);                                        // alpha0 = prevalence; u = pathway deviations
  }else{
    u_centered=rep_vector(0,P);                                  // No pathway: every gene shares the same prior log-odds alpha0
  }
  pi=inv_logit(alpha0+u_centered);                               // Convert log-odds to probability

  // === Spatial evidence for each gene ===
  if(use_shared_Phi==1){
    // Collapsed over w (exact): only r_eff terms per gene.
    matrix[r_eff,P] C=G'*z_ll;                                   // Project each gene's latent residuals onto the orthonormal basis
    for(j in 1:P){                                               // Loop over all P genes
      vector[r_eff] a=square(sig_eta_gs[j]/sigma_sd[j])*lam;     // Signal-to-noise ratio along each basis direction
      // a/(1+a) written as 1 - 1/(1+a): identical, but stays finite if a overflows in early warmup
      logLR[j]=-0.5*sum(log1p(a)-square(C[,j]).*(1-inv(1+a)));   // Log-likelihood ratio of SE vs non-SE for gene j
    }
  }else{
    for(j in 1:P){                                               // Loop over all P genes
      matrix[N,r] Phi_j;                                         // Gene-specific basis matrix
      for(i in 1:N){                                             // Loop over spots
        for(b in 1:r){                                           // Loop over knots
          Phi_j[i,b]=exp(-square(D[i,b])/(2*square(ell_gs[j])));  // RBF basis using gene j's own lengthscale
        }
      }
      vector[N] e_j=Phi_j*w[,j];                                 // Raw spatial effect for gene j
      vector[N] eta_j=sig_eta_gs[j]*(e_j-mean(e_j));             // centred over spots
      logLR[j]=normal_lpdf(z_ll[,j]|eta_j/sigma_sd[j],1)-std_normal_lpdf(z_ll[,j]);  // Log-likelihood ratio of SE vs non-SE for gene j
    }
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
  mu0~normal(mu_intercept,sd_intercept);                         // Prior on global intercept: mu0 ~ Normal(mu_intercept, sd_intercept)
  sigma_sd~normal(mu_noise,sd_noise);                            // Prior on residual SD for each gene; larger sd_noise allows more noise
  sig_eta_gs~normal(mu_gp_amplitude,sd_gp_amplitude);            // Prior on GP amplitude for each gene; controls expected strength of spatial signal
  to_vector(w)~std_normal();                                     // Prior on GP weights; standard normal on each knot weight (empty when use_shared_Phi = 1)
  ell_gs~lognormal(mu_gp_lengthscale,sd_gp_lengthscale);         // Prior on GP lengthscale for each gene
  sigma_beta~normal(mu_beta_sig,sd_beta_sig);                    // Prior on gene-level effect SD
  Beta~normal(0,sigma_beta);                                     // Prior on gene-level effects; independent across genes
  alpha0~normal(alpha0_mean,alpha0_sd);                          // Prior on global SE log-odds
  if(use_pathway==1){
    sigma_u[1]~normal(mu_u_sig,sd_u_sig);                        // Prior on pathway effect SD; controls how far genes can deviate from alpha0
  }
  if(use_pathway==1&&fix_rho_zero==0){
    rho_u_free[1]~beta(shape_rho_u,rate_rho_u);                  // Prior on CAR dependence; larger shape_rho_u pushes rho_u toward 1 (stronger borrowing)
  }

  // ================================================
  // SECTION 2: PATHWAY CAR PRIOR FOR u_j
  // ================================================
  // Genes in the same pathway are "neighbours" and share information on
  // their SE log-odds through the pathway effects u_j.
  // Leroux CAR (use_car_lx = 1):
  //   u_j | u_-j ~ N( rho_u sum_{l in N(j)} u_l / (1 - rho_u + rho_u n_j),
  //                   1 / (tau_u (1 - rho_u + rho_u n_j)) )
  // Sparse proper CAR (use_car_lx = 0):
  //   u_j | u_-j ~ N( rho_u sum_{l in N(j)} u_l / n_j,  1 / (tau_u n_j) )
  if(use_pathway==1){
    if(P_car>0){                                                 // Genes with at least one pathway neighbour
      if(use_car_lx==1){
        target+=leroux_car_lpdf(u[car_index]|sigma_u[1],rho_u,W_sparse,D_sparse,lap_eig,P_car,W_n);  // Leroux CAR log density
      }else{
        target+=sparse_car_lpdf(u[car_index]|sigma_u[1],rho_u,W_sparse,D_sparse,car_eig,P_car,W_n);  // Sparse proper CAR log density
      }
    }
    if(P_iso>0){
      u[iso_index]~normal(0,sigma_u[1]);                         // genes with no pathway neighbours
    }
  }

  // =============================================================
  // SECTION 3: LIKELIHOOD WITH Z_j MARGINALIZED
  // =============================================================
  // This is the core of the model. For each gene j:
  //   1. The observed counts Y[,j] follow a Poisson with rate N_i * lambda_ij
  //   2. loglambda[,j] = mu0 + Beta[j] + sigma_sd[j] * z_ll[,j]
  //   3. Non-SE: z ~ N(0, I);  SE: z ~ N(0, I + t_j^2 Phi_c Phi_c')
  //   4. log p(z_j) = log N(z_j | 0, I) + log_mix(pi_j, logLR_j, 0)
  for(j in 1:P){                                                 // Loop over all P genes
    vector[N] ll_j=mu0+Beta[j]+sigma_sd[j]*z_ll[,j];             // Latent log-rates for gene j across all spots
    target+=poisson_log_lpmf(Y[,j]|log_N_i+ll_j);                // Y[,j] ~ PoissonLog(log(N_i) + loglambda[,j])
    target+=std_normal_lpdf(z_ll[,j])+log_mix(pi[j],logLR[j],0); // Mixture over SE / non-SE with Z_j summed out
  }
}



// ---------------------------------------------------------------
// GENERATED QUANTITIES BLOCK - POSTERIOR CLASSIFICATION
// ---------------------------------------------------------------
// Computed AFTER each HMC draw, using the sampled parameter values.
// These are derived quantities for posterior inference.
// ---------------------------------------------------------------
generated quantities {
  // prob_SE[j] = P(Z_j=1 | rest) = inv_logit( alpha0 + u_centered[j] + logLR[j] )
  //            = pathway prior log-odds + spatial evidence
  array[P] int<lower=0,upper=1> Z;                               // Posterior gene-level classification: Z[j] = 0 non-SE; Z[j] = 1 SE
  vector<lower=0,upper=1>[P] prob_SE;                            // Posterior probability that gene j is SE
  for(j in 1:P){                                                 // Loop over all P genes
    prob_SE[j]=inv_logit(alpha0+u_centered[j]+logLR[j]);         // Combine pathway prior and spatial evidence
    Z[j]=bernoulli_rng(prob_SE[j]);                              // Draw Z[j]; averaging over draws gives P(gene j is SE)
  }
}
