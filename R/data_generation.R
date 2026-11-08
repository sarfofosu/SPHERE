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
#' @keywords internal
generate_spatial_spots <- function(num_spots,x1=0,x2=1,y1=0,y2=1,seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  data.frame(x=runif(num_spots,x1,x2),y=runif(num_spots,y1,y2))
}


## ------------------------------------------------------------
## Function to generate overlapping pathway membership
## ------------------------------------------------------------

#' Generate Overlapping Pathway Membership
#'
#' Generates a binary gene-by-pathway membership matrix in which genes may
#' belong to more than one pathway. A value of 1 means that the gene belongs
#' to that pathway. Every pathway is guaranteed to contain at least two genes.
#'
#' @param num_genes Integer. Number of genes (at least 2).
#' @param num_pathways Integer. Number of pathways (at least 1).
#' @param min_membership Integer. Minimum number of pathways per gene
#'   (default: 1).
#' @param max_membership Integer. Maximum number of pathways per gene; cannot
#'   exceed \code{num_pathways} (default: 3).
#' @param prob_multi Numeric in \eqn{[0, 1]} or NULL. Probability that a gene
#'   belongs to more than \code{min_membership} pathways. If NULL, the number
#'   of pathways per gene is drawn uniformly between \code{min_membership} and
#'   \code{max_membership} (default: NULL).
#' @param balanced Logical. If TRUE, each gene's first pathway is assigned in
#'   turn so pathway sizes are balanced (default: FALSE).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A binary matrix (\code{num_genes} x \code{num_pathways}) with row
#'   names \code{Gene1, Gene2, ...} and column names
#'   \code{Pathway1, Pathway2, ...}.
#'
#' @importFrom stats runif
#' @export
generate_pathway_membership <- function(num_genes,num_pathways,min_membership=1,max_membership=3,
                                        prob_multi=NULL,balanced=FALSE,seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  if(num_genes<2) stop("num_genes must be at least 2.")
  if(num_pathways<1) stop("num_pathways must be at least 1.")
  if(min_membership<1) stop("min_membership must be at least 1.")
  if(max_membership<min_membership) stop("max_membership must be >= min_membership.")
  if(max_membership>num_pathways) stop("max_membership cannot exceed num_pathways.")
  if(!is.null(prob_multi)&&(prob_multi<0||prob_multi>1)) stop("prob_multi must be in [0,1].")
  M <- matrix(0,num_genes,num_pathways)
  if(balanced) first_pw <- sample(rep_len(seq_len(num_pathways),num_genes))
  for(j in 1:num_genes){
    if(is.null(prob_multi)||max_membership==min_membership){
      k <- if(max_membership==min_membership) min_membership else sample(min_membership:max_membership,1)
    }else{
      # min_membership with prob 1-prob_multi, otherwise uniform on the larger sizes
      extra <- (min_membership+1):max_membership
      k <- if(runif(1)<prob_multi) extra[sample.int(length(extra),1)] else min_membership
    }
    if(balanced){
      rest <- setdiff(seq_len(num_pathways),first_pw[j])
      idx <- c(first_pw[j],rest[sample.int(length(rest),k-1)])
    }else{
      idx <- sample(seq_len(num_pathways),k)
    }
    M[j,idx] <- 1
  }
  for(k in 1:num_pathways){
    if(sum(M[,k])<2){
      idx <- sample(seq_len(num_genes),2)
      M[idx,k] <- 1
    }
  }
  rownames(M) <- paste0("Gene",seq_len(num_genes))
  colnames(M) <- paste0("Pathway",seq_len(num_pathways))
  M
}


## ------------------------------------------------------------
## Function to generate exclusive SE and non-SE pathways
## ------------------------------------------------------------

#' Generate Exclusive SE and Non-SE Pathway Membership
#'
#' Creates a controlled gene-pathway membership matrix together with the
#' true SE indicators, where selected pathways contain only SE genes and the
#' remaining pathways contain only non-SE genes. Used with
#' \code{se_design = "exclusive_pathways"} in \code{\link{make_sim_design}}.
#'
#' @param num_genes Integer. Number of genes.
#' @param num_pathways Integer. Number of pathways.
#' @param n_se Integer or NULL. Number of SE genes. If NULL,
#'   \code{round(prop_se * num_genes)} (default: NULL).
#' @param prop_se Numeric. Proportion of SE genes, used when \code{n_se} is
#'   NULL (default: 0.4).
#' @param n_se_pathways Integer or NULL. Number of SE-only pathways
#'   (default: NULL).
#' @param prop_se_pathways Numeric or NULL. Proportion of SE-only pathways,
#'   used when \code{n_se_pathways} is NULL (default: NULL).
#' @param se_pathways Integer or character vector, or NULL. Explicit SE-only
#'   pathways. Overrides \code{n_se_pathways} and \code{prop_se_pathways}
#'   (default: NULL).
#' @param min_membership Integer. Minimum number of pathways per gene
#'   (default: 1).
#' @param max_membership Integer. Maximum number of pathways per gene
#'   (default: 4).
#' @param balanced Logical. If TRUE, extra memberships go to the least-used
#'   pathways (default: TRUE).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A named list with elements \code{M} (membership matrix), \code{Z}
#'   (true SE indicators, 1 = SE), \code{se_idx}, \code{nonse_idx},
#'   \code{se_pathways}, \code{nonse_pathways}, \code{n_se}, \code{n_nonse},
#'   \code{n_se_pathways} and \code{prop_se_pathways}.
#'
#' @export
generate_exclusive_pathway_membership <- function(num_genes, num_pathways, n_se = NULL,
                                                  prop_se = 0.4, n_se_pathways = NULL, prop_se_pathways = NULL, se_pathways = NULL,
                                                  min_membership = 1, max_membership = 4, balanced = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- num_genes
  K <- num_pathways
  if (is.null(n_se)) n_se <- round(prop_se * P)
  n_se <- as.integer(n_se)
  if (is.null(se_pathways)) {
    if (is.null(n_se_pathways)) {
      if (is.null(prop_se_pathways)) stop("Specify n_se_pathways, prop_se_pathways or se_pathways.")
      n_se_pathways <- round(prop_se_pathways * K)
    }
    n_se_pathways <- as.integer(n_se_pathways)
    if (length(n_se_pathways) != 1 || is.na(n_se_pathways) || n_se_pathways < 1 || n_se_pathways >= K) stop("n_se_pathways must be between 1 and num_pathways-1.")
    se_pathways <- sort(sample.int(K, n_se_pathways))
  } else {
    if (is.character(se_pathways)) se_pathways <- match(se_pathways, paste0("Pathway", seq_len(K)))
    se_pathways <- sort(unique(as.integer(se_pathways)))
  }
  if (length(se_pathways) == 0 || anyNA(se_pathways) || any(se_pathways < 1 | se_pathways > K) || length(se_pathways) >= K) stop("Invalid se_pathways.")
  if (length(n_se) != 1 || is.na(n_se) || n_se < 1 || n_se >= P) stop("n_se must be between 1 and num_genes-1.")
  if (length(min_membership) != 1 || length(max_membership) != 1 || is.na(min_membership) || is.na(max_membership) || min_membership < 1 || max_membership < min_membership || min_membership != as.integer(min_membership) || max_membership != as.integer(max_membership)) stop("Invalid membership limits.")
  nonse_pathways <- setdiff(seq_len(K), se_pathways)
  n_nonse <- P - n_se
  if (n_se < length(se_pathways) || n_nonse < length(nonse_pathways)) stop("Each pathway needs at least one gene. Increase the group size or reduce its number of pathways.")
  if (min_membership > length(se_pathways) || min_membership > length(nonse_pathways)) stop("min_membership exceeds the number of pathways available to one gene group.")

  # Safe sampling: sample(x,n) with a single number x draws from 1:x, so pick from the vector itself.
  pick <- function(x, n) x[sample.int(length(x), n)]

  # Allocate genes to pathways within their own SE class.
  make_group <- function(n_genes, pathways) {
    k <- length(pathways)
    A <- matrix(0L, n_genes, k)

    # Guarantee that every pathway contains at least one gene.
    order_genes <- sample.int(n_genes)
    order_pathways <- sample.int(k)
    for (i in seq_len(k)) A[order_genes[i], order_pathways[i]] <- 1L

    # Ensure every gene has at least min_membership pathways.
    for (i in seq_len(n_genes)) {
      need <- min_membership - sum(A[i, ])
      if (need > 0) {
        available <- which(A[i, ] == 0)
        selected <- pick(available, need)
        A[i, selected] <- 1L
      }
    }

    # Assign additional memberships up to max_membership.
    upper <- min(max_membership, k)
    if (upper > min_membership) {
      for (i in seq_len(n_genes)) {
        target <- sample(min_membership:upper, 1)
        need <- target - sum(A[i, ])
        if (need > 0) {
          available <- which(A[i, ] == 0)
          if (balanced) {
            selected <- integer(0)
            for (h in seq_len(need)) {
              remaining <- setdiff(available, selected)
              loads <- colSums(A)[remaining]
              candidates <- remaining[loads == min(loads)]
              selected <- c(selected, pick(candidates, 1))
            }
          } else {
            selected <- pick(available, need)
          }
          A[i, selected] <- 1L
        }
      }
    }
    colnames(A) <- paste0("Pathway", pathways)
    A
  }

  # Randomly assign exactly n_se genes as SE.
  Z <- integer(P)
  Z[sample.int(P, n_se)] <- 1L
  se_idx <- which(Z == 1)
  nonse_idx <- which(Z == 0)
  M <- matrix(0L, P, K)
  colnames(M) <- paste0("Pathway", seq_len(K))
  rownames(M) <- paste0("Gene", seq_len(P))
  M[se_idx, se_pathways] <- make_group(n_se, se_pathways)
  M[nonse_idx, nonse_pathways] <- make_group(n_nonse, nonse_pathways)
  list(M = M, Z = Z, se_idx = se_idx, nonse_idx = nonse_idx, se_pathways = se_pathways, nonse_pathways = nonse_pathways, n_se = n_se, n_nonse = n_nonse, n_se_pathways = length(se_pathways), prop_se_pathways = length(se_pathways) / K)
}


## ------------------------------------------------------------
## Function to generate exclusive pathways with non-SE genes mixed in
## ------------------------------------------------------------

#' Generate Exclupath Pathway Membership
#'
#' Creates an "exclusive pathways with non-SE genes mixed in" design. A set of
#' pathways are SE pathways and the rest are non-SE pathways. SE genes appear
#' only in SE pathways. Each SE pathway holds exactly \code{se_per_pathway} SE
#' genes plus \code{nonse_per_se_pathway} non-SE genes, so SE pathways are
#' never SE-only. Non-SE pathways hold only non-SE genes with balanced sizes.
#' Used with \code{se_design = "exclupath"} in \code{\link{make_sim_design}}.
#'
#' @param num_genes Integer. Number of genes.
#' @param num_pathways Integer. Number of pathways.
#' @param n_se Integer. Number of SE genes.
#' @param n_se_pathways Integer or NULL. Number of SE pathways (default: NULL).
#' @param prop_se_pathways Numeric. Proportion of SE pathways, used when
#'   \code{n_se_pathways} is NULL (default: 0.5).
#' @param se_pathways Integer vector or NULL. Explicit SE pathways
#'   (default: NULL).
#' @param se_per_pathway Integer or NULL. SE genes in each SE pathway. Must
#'   satisfy \code{n_se <= se_per_pathway * n_se_pathways}. If NULL,
#'   \code{ceiling(n_se / n_se_pathways)} (default: NULL).
#' @param nonse_per_se_pathway Integer (at least 1). Non-SE genes added to
#'   each SE pathway (default: 1).
#' @param contam_bridge Logical. If FALSE, the non-SE genes placed in SE
#'   pathways belong only to that SE pathway; if TRUE they also keep a non-SE
#'   pathway, linking SE and non-SE pathways (default: FALSE).
#' @param nonse_prob_multi Numeric. Probability that a non-SE gene also joins
#'   a second non-SE pathway (default: 0).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A named list with elements \code{M} (membership matrix), \code{Z}
#'   (true SE indicators, 1 = SE), \code{se_idx}, \code{nonse_idx},
#'   \code{se_pathways}, \code{nonse_pathways}, \code{contaminants},
#'   \code{se_per_pathway}, \code{nonse_per_se_pathway}, \code{n_se} and
#'   \code{n_se_pathways}.
#'
#' @importFrom stats runif setNames
#' @export
generate_exclupath_membership <- function(num_genes, num_pathways, n_se, n_se_pathways = NULL, prop_se_pathways = 0.5,
                                          se_pathways = NULL, se_per_pathway = NULL, nonse_per_se_pathway = 1,
                                          contam_bridge = FALSE, nonse_prob_multi = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- as.integer(num_genes); K <- as.integer(num_pathways); n_se <- as.integer(n_se)
  pick <- function(x, n) x[sample.int(length(x), n)]
  # which pathways are SE pathways
  if (is.null(se_pathways)) {
    if (is.null(n_se_pathways)) n_se_pathways <- round(prop_se_pathways * K)
    se_pathways <- sort(sample.int(K, n_se_pathways))
  } else se_pathways <- sort(unique(as.integer(se_pathways)))
  k_se <- length(se_pathways); nonse_pathways <- setdiff(seq_len(K), se_pathways); k_non <- length(nonse_pathways)
  # [CHANGED 3] non-SE pathways are no longer required (all pathways can be SE pathways), and
  # every SE pathway must hold at least one non-SE gene.
  if (k_se < 1) stop("Need at least one SE pathway.", call. = FALSE)
  if (nonse_per_se_pathway < 1 || nonse_per_se_pathway != round(nonse_per_se_pathway))
    stop("nonse_per_se_pathway must be a whole number >= 1 (every SE pathway needs a non-SE gene).", call. = FALSE)
  if (is.null(se_per_pathway)) se_per_pathway <- ceiling(n_se / k_se)
  if (se_per_pathway * k_se < n_se) stop("se_per_pathway x n_se_pathways must be >= n_se so every SE gene has a pathway.", call. = FALSE)
  if (se_per_pathway > n_se) stop("se_per_pathway cannot exceed n_se.", call. = FALSE)
  n_non <- P - n_se
  n_contam <- nonse_per_se_pathway * k_se
  if ((if (contam_bridge) 0 else n_contam) + k_non > n_non) stop("Not enough non-SE genes for the SE pathways and non-SE pathways.", call. = FALSE)
  # [CHANGED 4] with no non-SE pathways, every non-SE gene must go inside the SE pathways
  if (k_non == 0 && (contam_bridge || n_contam != n_non))
    stop(sprintf("All %d pathways are SE pathways, so the %d non-SE genes must all go into them: set nonse_per_se_pathway = %g.",
                 K, n_non, n_non / k_se), call. = FALSE)

  Z <- integer(P); Z[sample.int(P, n_se)] <- 1L
  se_idx <- which(Z == 1); non_idx <- which(Z == 0)
  M <- matrix(0L, P, K, dimnames = list(paste0("Gene", seq_len(P)), paste0("Pathway", seq_len(K))))

  # SE genes -> SE pathways: exactly se_per_pathway each, spreading genes evenly
  # (always pick the least-used SE genes not yet in this pathway).
  use <- setNames(integer(n_se), se_idx)
  for (k in se_pathways[sample.int(k_se)]) {
    cand <- se_idx
    chosen <- integer(0)
    for (h in seq_len(se_per_pathway)) {
      left <- setdiff(cand, chosen)
      u <- use[as.character(left)]
      g <- pick(left[u == min(u)], 1)
      chosen <- c(chosen, g); use[as.character(g)] <- use[as.character(g)] + 1L
    }
    M[chosen, k] <- 1L
  }

  # Non-SE genes placed in SE pathways ("contamination"), the same number in each.
  contam <- if (n_contam > 0) pick(non_idx, n_contam) else integer(0)
  if (n_contam > 0) M[cbind(contam, rep(se_pathways, each = nonse_per_se_pathway))] <- 1L

  # Remaining non-SE genes -> non-SE pathways, balanced (dealt out in turn).
  rest <- if (contam_bridge) non_idx else setdiff(non_idx, contam)
  rest <- rest[sample.int(length(rest))]
  M[cbind(rest, rep_len(nonse_pathways[sample.int(k_non)], length(rest)))] <- 1L
  # Optional second non-SE pathway for some non-SE genes (to the least-loaded one).
  if (nonse_prob_multi > 0 && k_non > 1) for (g in rest) if (runif(1) < nonse_prob_multi) {
    free <- nonse_pathways[M[g, nonse_pathways] == 0]
    load <- colSums(M)[free]; M[g, pick(free[load == min(load)], 1)] <- 1L
  }
  list(M = M, Z = Z, se_idx = se_idx, nonse_idx = non_idx, se_pathways = se_pathways, nonse_pathways = nonse_pathways,
       contaminants = contam, se_per_pathway = se_per_pathway, nonse_per_se_pathway = nonse_per_se_pathway,
       n_se = n_se, n_se_pathways = k_se)
}


## ------------------------------------------------------------
## Function to simulate one Gaussian Process spatial pattern
## ------------------------------------------------------------

#' Simulate One Gaussian Process Spatial Pattern
#'
#' Draws a single spatial effect from a zero-mean Gaussian Process with a
#' squared-exponential (RBF) kernel
#' \eqn{K = \sigma_\eta^2 \exp(-d^2 / (2 \ell^2))}.
#'
#' @param spots A numeric matrix or data frame of spatial coordinates.
#' @param sigma_eta Numeric. GP amplitude; controls the strength of the
#'   spatial signal.
#' @param ell Numeric. GP lengthscale; small values give localized patterns,
#'   large values give smoother patterns.
#' @param nugget Numeric. Nugget added to the covariance diagonal for
#'   numerical stability (default: 1e-8).
#'
#' @return A numeric vector of length \code{nrow(spots)}.
#'
#' @importFrom stats dist rnorm
#' @keywords internal
simulate_full_gp <- function(spots,sigma_eta,ell,nugget=1e-8){
  n <- nrow(spots)
  if(sigma_eta<=0) stop("sigma_eta must be positive.")
  if(ell<=0) stop("ell must be positive.")
  d2 <- as.matrix(dist(spots))^2
  K <- sigma_eta^2*exp(-d2/(2*ell^2))
  diag(K) <- diag(K)+nugget
  R <- chol(K)
  eta <- as.vector(t(R)%*%rnorm(n))
  eta
}


## ------------------------------------------------------------
## Function to simulate pathway effects u (CAR)
## ------------------------------------------------------------

#' Simulate Pathway Effects from a CAR Prior
#'
#' Simulates gene-level pathway effects \eqn{u} on the SE log-odds from a
#' Conditional Autoregressive (CAR) prior defined on the pathway network.
#' With \code{car_type = "leroux"},
#' \eqn{u \sim N(0, \sigma_u^2 [(1-\rho_u) I + \rho_u (D - W)]^{-1})};
#' with \code{car_type = "proper"}, the precision is \eqn{D - \rho_u W}.
#'
#' @param num_genes Integer. Number of genes.
#' @param gene_grp Integer vector or NULL. Single pathway label per gene.
#'   Supply this or \code{pathway_membership}, not both (default: NULL).
#' @param pathway_membership Binary matrix or NULL. Gene-by-pathway
#'   membership (default: NULL).
#' @param rho_u Numeric in \eqn{[0, 1)}. CAR dependence (default: 0.8).
#' @param sigma_u Numeric. Scale of the pathway effects (default: 1).
#' @param car_type Character. \code{"leroux"} (default) or \code{"proper"}.
#' @param standardize Logical. If TRUE, \eqn{u} is scaled so its average
#'   centred marginal SD equals \code{sigma_u} (default: TRUE).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A named list including \code{u}, the network \code{W}, degree
#'   matrix \code{D}, base precision \code{Q0}, covariance \code{Sigma_u},
#'   \code{scale}, \code{degree}, \code{membership}, \code{shared_pathways}
#'   and the settings used.
#'
#' @importFrom stats rnorm
#' @keywords internal
create_car_u<-function(num_genes,gene_grp=NULL,pathway_membership=NULL,rho_u=0.8,sigma_u=1,
                       car_type=c("leroux","proper"),standardize=TRUE,seed=NULL){
  car_type<-match.arg(car_type)
  if(!is.null(seed))set.seed(seed)
  if(rho_u<0||rho_u>=1)stop("rho_u must satisfy 0 <= rho_u < 1.")
  if(sigma_u<=0)stop("sigma_u must be positive.")
  net<-build_pathway_network(num_genes,gene_grp,pathway_membership)
  W<-net$W
  degree<-net$degree
  D<-diag(degree,nrow=num_genes)
  if(car_type=="leroux"){
    Q0<-(1-rho_u)*diag(num_genes)+rho_u*(D-W)
  }else{
    Q0<-D-rho_u*W
  }
  iso<-which(degree==0)
  if(length(iso)>0)diag(Q0)[iso]<-1
  Q0<-(Q0+t(Q0))/2
  R<-tryCatch(chol(Q0),error=function(e)stop("CAR precision is not positive definite.",call.=FALSE))
  Sigma0<-chol2inv(R)
  if(standardize){
    C<-diag(num_genes)-matrix(1/num_genes,num_genes,num_genes)
    scale<-sigma_u/sqrt(mean(diag(C%*%Sigma0%*%C)))
  }else{
    scale<-sigma_u
  }
  u<-scale*as.numeric(backsolve(R,rnorm(num_genes)))
  Sigma_u<-scale^2*Sigma0
  list(u=u,W=W,D=D,Q0=Q0,Sigma_u=Sigma_u,scale=scale,degree=degree,membership=net$membership,
       shared_pathways=net$shared_pathways,rho_u=rho_u,sigma_u=sigma_u,car_type=car_type,standardize=standardize)
}


## ------------------------------------------------------------
## Function to assign SE status Z from pathway effects
## ------------------------------------------------------------

#' Assign Pathway-Informed SE Status
#'
#' Generates SE status through the chain: pathway network -> pathway effect
#' \eqn{u} (CAR) -> SE probability \eqn{\pi_j = \mathrm{logit}^{-1}(\alpha_0 + u_j)}
#' -> \eqn{Z_j \sim \mathrm{Bernoulli}(\pi_j)}. \eqn{\alpha_0} is chosen so the
#' average SE probability equals \code{base_se_prob}. Also returns diagnostics
#' of how much detection information the pathway carries.
#'
#' @inheritParams create_car_u
#' @param base_se_prob Numeric in \eqn{(0, 1)}. Target average SE probability
#'   (default: 0.2).
#'
#' @return A named list with \code{Z} (1 = SE), \code{se_idx},
#'   \code{nonse_idx}, \code{n_se}, \code{alpha0}, \code{u},
#'   \code{u_centered}, \code{pi}, the network quantities from
#'   \code{create_car_u}, \code{nbr_corr} (correlation between \eqn{u_j} and
#'   its neighbours' mean) and \code{sd_logit_pi} (SD of the SE log-odds).
#'
#' @importFrom stats plogis uniroot rbinom cor sd qlogis
#' @keywords internal
assign_pathway_se_status <- function(num_genes, gene_grp = NULL, pathway_membership = NULL,
                                     base_se_prob = 0.2, rho_u = 0.8, sigma_u = 1,
                                     car_type = c("leroux", "proper"), standardize = TRUE, seed = NULL) {
  car_type <- match.arg(car_type)
  if (!is.null(seed)) set.seed(seed)
  if (base_se_prob <= 0 || base_se_prob >= 1) stop("base_se_prob must satisfy 0 < base_se_prob < 1.")
  # alpha0 <- qlogis(base_se_prob)
  # car <- create_car_u(num_genes, gene_grp, pathway_membership, rho_u, sigma_u, car_type = car_type, standardize = standardize)
  # u <- car$u
  # u_centered <- u - mean(u)
  # pi <- plogis(alpha0 + u_centered)
  car <- create_car_u(num_genes, gene_grp, pathway_membership, rho_u, sigma_u, car_type = car_type, standardize = standardize)
  u <- car$u
  u_centered <- u - mean(u)
  # Choose alpha0 so the average SE probability equals base_se_prob (whatever sigma_u and rho_u are)
  alpha0 <- uniroot(function(a) mean(plogis(a + u_centered)) - base_se_prob, c(-30, 30))$root
  pi <- plogis(alpha0 + u_centered)
  Z <- rbinom(num_genes, 1, pi)
  se_idx <- which(Z == 1)
  nonse_idx <- which(Z == 0)
  has_nb <- car$degree > 0
  nbr_mean <- rep(NA_real_, num_genes)
  nbr_mean[has_nb] <- as.vector(car$W %*% u_centered)[has_nb] / car$degree[has_nb]
  nbr_corr <- if (sum(has_nb) > 2) cor(u_centered[has_nb], nbr_mean[has_nb]) else NA_real_
  list(
    Z = Z, se_idx = se_idx, nonse_idx = nonse_idx, n_se = length(se_idx), alpha0 = alpha0, u = u, u_centered = u_centered, pi = pi,
    W = car$W, D = car$D, Q0 = car$Q0, Sigma_u = car$Sigma_u, degree = car$degree, membership = car$membership,
    shared_pathways = car$shared_pathways, rho_u = rho_u, sigma_u = sigma_u, car_type = car_type,
    nbr_corr = nbr_corr, sd_logit_pi = sd(qlogis(pi))
  )
}


## ------------------------------------------------------------
## Function to assign pattern strength and width for SE genes
## ------------------------------------------------------------

#' Assign GP Strength and Lengthscale for SE Genes
#'
#' Draws a GP amplitude and lengthscale for every SE gene. Non-SE genes get
#' amplitude 0 and lengthscale NA.
#'
#' @param P Integer. Number of genes.
#' @param se_idx Integer vector. Indices of SE genes.
#' @param signal_scenario Character. One of \code{"weak"},
#'   \code{"moderate"}, \code{"strong"}, \code{"mixed"} or \code{"custom"}
#'   (default: \code{"mixed"}).
#' @param sigma_eta_range Numeric vector of length 2. Amplitude range for
#'   \code{"custom"} (default: \code{c(0.4, 0.7)}).
#' @param weak_range Numeric vector of length 2. Amplitude range for weak
#'   signals (default: \code{c(0.20, 0.35)}).
#' @param moderate_range Numeric vector of length 2. Amplitude range for
#'   moderate signals (default: \code{c(0.35, 0.55)}).
#' @param strong_range Numeric vector of length 2. Amplitude range for
#'   strong signals (default: \code{c(0.55, 0.80)}).
#' @param ell_range Numeric vector of length 2. Lengthscale range
#'   (default: \code{c(0.15, 0.35)}).
#' @param signal_probs Named numeric vector (weak, moderate, strong) summing
#'   to 1. Class probabilities for \code{"mixed"}
#'   (default: \code{c(weak = 0.3, moderate = 0.3, strong = 0.4)}).
#'
#' @return A named list with \code{sigma_eta}, \code{ell} and
#'   \code{signal_class} (each of length \code{P}).
#'
#' @importFrom stats runif
#' @keywords internal
assign_gp_parameters <- function(P, se_idx, signal_scenario = "mixed", sigma_eta_range = c(0.4, 0.7),
                                 weak_range = c(0.20, 0.35), moderate_range = c(0.35, 0.55),
                                 strong_range = c(0.55, 0.80), ell_range = c(0.15, 0.35),
                                 signal_probs = c(weak = 0.3, moderate = 0.3, strong = 0.4)) {
  valid <- c("weak", "moderate", "strong", "mixed", "custom")
  if (!signal_scenario %in% valid) stop("signal_scenario must be weak, moderate, strong, mixed, or custom.")
  if (length(ell_range) != 2 || ell_range[1] <= 0 || ell_range[2] < ell_range[1]) stop("Invalid ell_range.")
  sigma_eta <- rep(0, P)
  ell <- rep(NA_real_, P)
  signal_class <- rep("NonSE", P)
  n_se <- length(se_idx)
  if (n_se == 0) {
    return(list(sigma_eta = sigma_eta, ell = ell, signal_class = signal_class))
  }
  if (signal_scenario == "weak") {
    sigma_eta[se_idx] <- runif(n_se, weak_range[1], weak_range[2])
    signal_class[se_idx] <- "weak"
  }
  if (signal_scenario == "moderate") {
    sigma_eta[se_idx] <- runif(n_se, moderate_range[1], moderate_range[2])
    signal_class[se_idx] <- "moderate"
  }
  if (signal_scenario == "strong") {
    sigma_eta[se_idx] <- runif(n_se, strong_range[1], strong_range[2])
    signal_class[se_idx] <- "strong"
  }
  if (signal_scenario == "custom") {
    sigma_eta[se_idx] <- runif(n_se, sigma_eta_range[1], sigma_eta_range[2])
    signal_class[se_idx] <- "custom"
  }
  if (signal_scenario == "mixed") {
    if (length(signal_probs) != 3 || any(signal_probs < 0) || abs(sum(signal_probs) - 1) > 1e-8) stop("signal_probs must contain weak, moderate, strong probabilities summing to 1.")
    cls <- sample(c("weak", "moderate", "strong"), n_se, replace = TRUE, prob = signal_probs)
    signal_class[se_idx] <- cls
    weak_idx <- se_idx[cls == "weak"]
    moderate_idx <- se_idx[cls == "moderate"]
    strong_idx <- se_idx[cls == "strong"]
    if (length(weak_idx) > 0) sigma_eta[weak_idx] <- runif(length(weak_idx), weak_range[1], weak_range[2])
    if (length(moderate_idx) > 0) sigma_eta[moderate_idx] <- runif(length(moderate_idx), moderate_range[1], moderate_range[2])
    if (length(strong_idx) > 0) sigma_eta[strong_idx] <- runif(length(strong_idx), strong_range[1], strong_range[2])
  }
  ell[se_idx] <- runif(n_se, ell_range[1], ell_range[2])
  list(sigma_eta = sigma_eta, ell = ell, signal_class = signal_class)
}


## ------------------------------------------------------------
## Function -- SPHERE model data generator
## ------------------------------------------------------------

#' Generate Model-Based Spatial Transcriptomics Count Data
#'
#' Simulates one spatial transcriptomics dataset consistent with the SPHERE
#' model. SE status is pathway-informed: pathway network -> pathway effect
#' \eqn{u} (CAR) -> \eqn{\pi_j} -> \eqn{Z_j \sim \mathrm{Bernoulli}(\pi_j)}. SE
#' genes receive a GP spatial pattern of random strength; non-SE genes get
#' none. Counts follow
#' \eqn{Y_{ij} \sim \mathrm{Poisson}(N_i \lambda_{ij})} with
#' \eqn{\log \lambda_{ij} = \mu_0 + \beta_j + \varepsilon_{ij} + \eta_{ij}}.
#'
#' @param num_spots Integer. Number of spatial spots.
#' @param num_genes Integer. Number of genes.
#' @param gene_grp Integer vector or NULL. Single pathway label per gene.
#'   Supply this or \code{pathway_membership}, not both (default: NULL).
#' @param pathway_membership Binary matrix or NULL. Gene-by-pathway
#'   membership, e.g. from \code{\link{generate_pathway_membership}}
#'   (default: NULL).
#' @param base_se_prob Numeric in \eqn{(0, 1)}. Target average SE probability
#'   (default: 0.2).
#' @param rho_u Numeric in \eqn{[0, 1)}. CAR dependence of the pathway
#'   effects (default: 0.8).
#' @param sigma_u Numeric. Scale of the pathway effects (default: 1).
#' @param sigma_beta Numeric. SD of the gene-level effects \eqn{\beta_j}
#'   (default: 0.3).
#' @param eps_sd Numeric. SD of the iid log-scale noise (default: 0.25).
#' @param signal_scenario Character. One of \code{"weak"},
#'   \code{"moderate"}, \code{"strong"}, \code{"mixed"} or \code{"custom"}
#'   (default: \code{"mixed"}).
#' @param weak_range Numeric vector of length 2. GP amplitude range for weak
#'   signals (default: \code{c(0.10, 0.25)}).
#' @param moderate_range Numeric vector of length 2. GP amplitude range for
#'   moderate signals (default: \code{c(0.25, 0.45)}).
#' @param strong_range Numeric vector of length 2. GP amplitude range for
#'   strong signals (default: \code{c(0.45, 0.70)}).
#' @param custom_range Numeric vector of length 2. GP amplitude range for
#'   \code{"custom"} (default: \code{c(0.4, 0.7)}).
#' @param ell_range Numeric vector of length 2. GP lengthscale range on the
#'   unit square; rescaled by the domain size (default: \code{c(0.15, 0.35)}).
#' @param signal_probs Named numeric vector (weak, moderate, strong) summing
#'   to 1. Class probabilities for \code{"mixed"}
#'   (default: \code{c(weak = 0.4, moderate = 0.3, strong = 0.3)}).
#' @param depth_median Numeric. Median library size (default: 5000).
#' @param depth_log_sd Numeric. SD of log library size (default: 0.35).
#' @param baseline_mean Numeric. Baseline expected count; sets
#'   \code{mu0 = log(baseline_mean / depth_median)} (default: 5).
#' @param x1 Numeric. Minimum x-coordinate (default: 0).
#' @param x2 Numeric. Maximum x-coordinate (default: 1).
#' @param y1 Numeric. Minimum y-coordinate (default: 0).
#' @param y2 Numeric. Maximum y-coordinate (default: 1).
#' @param nugget Numeric. Nugget added to the GP covariance diagonal
#'   (default: 1e-8).
#' @param car_type Character. \code{"leroux"} (default) or \code{"proper"}.
#' @param standardize_u Logical. Standardize the pathway effects so their
#'   average marginal SD equals \code{sigma_u} (default: TRUE).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#' @param verbose Logical. If TRUE, prints a simulation summary (default: TRUE).
#'
#' @return A named list including \code{Y} (count matrix), \code{spots},
#'   \code{N_i}, \code{mu0}, \code{Beta}, \code{epsilon}, \code{eta},
#'   \code{loglambda}, \code{lambda}, \code{Y_mean}, \code{Z} (1 = SE),
#'   \code{se_idx}, \code{nonse_idx}, \code{alpha0}, \code{u},
#'   \code{u_centered}, \code{pi}, the pathway network (\code{W}, \code{D},
#'   \code{Q0}, \code{Sigma_u}, \code{degree}, \code{pathway_membership},
#'   \code{shared_pathways}), \code{sigma_eta}, \code{ell},
#'   \code{signal_class}, diagnostics (\code{nbr_corr}, \code{sd_logit_pi})
#'   and \code{param}, a list of all simulation settings.
#'
#' @importFrom stats rnorm rpois median
#' @export
generate_genedata_model <- function(num_spots, num_genes, gene_grp = NULL, pathway_membership = NULL,
                                    base_se_prob = 0.2, rho_u = 0.8, sigma_u = 1, sigma_beta = 0.3,
                                    eps_sd = 0.25, signal_scenario = "mixed", weak_range = c(0.10, 0.25),
                                    moderate_range = c(0.25, 0.45), strong_range = c(0.45, 0.70),
                                    custom_range = c(0.4, 0.7), ell_range = c(0.15, 0.35),
                                    signal_probs = c(weak = 0.4, moderate = 0.3, strong = 0.3),
                                    depth_median = 5000, depth_log_sd = 0.35, baseline_mean = 5,
                                    x1 = 0, x2 = 1, y1 = 0, y2 = 1, nugget = 1e-8, car_type = c("leroux", "proper"),
                                    standardize_u = TRUE, seed = NULL, verbose = TRUE) {
  car_type <- match.arg(car_type)
  if (!is.null(seed)) set.seed(seed)
  n <- num_spots
  P <- num_genes
  spots <- generate_spatial_spots(n, x1, x2, y1, y2)
  coords <- as.matrix(spots[, c("x", "y")])
  domain_scale <- sqrt((x2 - x1) * (y2 - y1))
  ell_range_scaled <- ell_range * domain_scale
  N_i <- round(exp(rnorm(n, log(depth_median), depth_log_sd)))
  N_i <- pmax(N_i, 1L)

  # Generate pathway-informed SE status.
  status <- assign_pathway_se_status(P, gene_grp, pathway_membership, base_se_prob, rho_u, sigma_u,
                                     car_type = car_type, standardize = standardize_u
  )
  Z <- status$Z
  se_idx <- status$se_idx
  nonse_idx <- status$nonse_idx
  n_se <- status$n_se
  alpha0 <- status$alpha0
  u <- status$u
  u_centered <- status$u_centered
  pi <- status$pi

  # Generate independent gene-level baseline effects.
  Beta <- rnorm(P, 0, sigma_beta)

  gp <- assign_gp_parameters(P = P, se_idx = se_idx, signal_scenario = signal_scenario,
                             sigma_eta_range = custom_range, weak_range = weak_range,
                             moderate_range = moderate_range, strong_range = strong_range,
                             signal_probs = signal_probs, ell_range = ell_range_scaled)
  sigma_eta <- gp$sigma_eta
  ell <- gp$ell
  signal_class <- gp$signal_class
  epsilon <- matrix(rnorm(n * P, 0, eps_sd), n, P)
  eta <- matrix(0, n, P)
  if (n_se > 0) {
    for (j in se_idx) eta[, j] <- simulate_full_gp(coords, sigma_eta[j], ell[j], nugget)
  }
  mu0 <- log(baseline_mean / depth_median)
  loglambda <- matrix(mu0, n, P) + matrix(Beta, n, P, byrow = TRUE) + epsilon + eta
  lambda <- exp(loglambda)
  Y_mean <- lambda * N_i
  Y <- matrix(rpois(n * P, as.vector(Y_mean)), n, P)
  gene_names <- character(P)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx] <- paste0("Gene", se_idx, "_SE_", signal_class[se_idx])
  colnames(Y) <- gene_names
  colnames(loglambda) <- gene_names
  colnames(eta) <- gene_names
  if (verbose) {
    cat("---- Pathway-Informed SPHERE Simulation ----\n")
    cat("Spots:", n, " Genes:", P, "\n")
    cat("SE genes:", n_se, "(", round(100 * n_se / P, 1), "%)\n")
    cat("Pathways:", ncol(status$membership), "\n")
    cat("Mean pathways per gene:", round(mean(rowSums(status$membership)), 2), "\n")
    cat("Genes in multiple pathways:", sum(rowSums(status$membership) > 1), "\n")
    cat("Mean network degree:", round(mean(status$degree), 2), "\n")
    cat("CAR type:", car_type, " standardized u:", standardize_u, "\n")
    cat("rho_u:", rho_u, " sigma_u:", sigma_u, "\n")
    cat("alpha0:", round(alpha0, 4), " baseline SE probability:", base_se_prob, "\n")
    cat("Mean pi:", round(mean(pi), 4), "\n")
    cat("Range pi:", paste(round(range(pi), 4), collapse = "-"), "\n")
    cat("SD of SE log-odds across genes:", round(status$sd_logit_pi, 3), "\n")
    cat(
      "corr(u_j, neighbour mean u):", round(status$nbr_corr, 3),
      " <- pathway information available to the CAR prior\n"
    )
    cat("Realized SE proportion:", round(mean(Z), 4), "\n")
    cat("sigma_beta:", sigma_beta, "\n")
    cat("Signal scenario:", signal_scenario, "\n")
    cat("Scaled ell range:", paste(round(ell_range_scaled, 3), collapse = "-"), "\n")
    if (n_se > 0) cat("Signal classes:", paste(names(table(signal_class[se_idx])), table(signal_class[se_idx]), collapse = ", "), "\n")
    cat("eps_sd:", eps_sd, "\n")
    cat("Depth median:", median(N_i), "\n")
    cat("Mean count:", round(mean(Y), 3), "\n")
  }
  list(
    Y = Y, spots = spots, N_i = N_i, mu0 = mu0, Beta = Beta, epsilon = epsilon,
    eta = eta, loglambda = loglambda, lambda = lambda, Y_mean = Y_mean, Z = Z,
    se_idx = se_idx, nonse_idx = nonse_idx, alpha0 = alpha0, u = u,
    u_centered = u_centered, pi = pi, W = status$W, D = status$D,
    Q0 = status$Q0, Sigma_u = status$Sigma_u, degree = status$degree,
    pathway_membership = status$membership, shared_pathways = status$shared_pathways,
    sigma_eta = sigma_eta, ell = ell, signal_class = signal_class, domain_scale = domain_scale,
    ell_range_scaled = ell_range_scaled,
    nbr_corr = status$nbr_corr, sd_logit_pi = status$sd_logit_pi,
    param = list(n = n, p = P, base_se_prob = base_se_prob, rho_u = rho_u, sigma_u = sigma_u,
                 sigma_beta = sigma_beta, eps_sd = eps_sd, signal_scenario = signal_scenario,
                 weak_range = weak_range, moderate_range = moderate_range, strong_range = strong_range,
                 signal_probs = signal_probs, ell_range = ell_range, domain_scale = domain_scale,
                 ell_range_scaled = ell_range_scaled, car_type = car_type, standardize_u = standardize_u, seed = seed)
  )
}


## ------------------------------------------------------------
## Function to turn n_hot into a number of pathways
## ------------------------------------------------------------

#' Resolve the Number of Hot Pathways
#'
#' Turns \code{n_hot} (\code{"all"}, \code{"half"} or a number) into a number
#' of pathways.
#'
#' @param n_hot Character (\code{"all"} or \code{"half"}) or integer.
#' @param K Integer. Total number of pathways.
#'
#' @return An integer between 1 and \code{K}.
#'
#' @keywords internal
resolve_n_hot <- function(n_hot, K) {
  if (is.character(n_hot)) {
    n_hot <- switch(match.arg(n_hot, c("all", "half")),
                    "all" = K,
                    "half" = ceiling(K / 2)
    )
  }
  n_hot <- as.integer(n_hot)
  if (length(n_hot) != 1 || is.na(n_hot) || n_hot < 1 || n_hot > K) {
    stop("n_hot must be 'all', 'half' or a number between 1 and ", K, " (the number of pathways).", call. = FALSE)
  }
  n_hot
}


## ------------------------------------------------------------
## Function to place SE genes in hot pathways
## ------------------------------------------------------------

#' Place SE Genes Inside Hot Pathways
#'
#' Places \code{n_se} SE genes inside the hot pathways. Genes that belong
#' only to hot pathways are used first. If there are not enough, then with
#' \code{strict = TRUE} the function stops; with \code{strict = FALSE} it also
#' uses hot-pathway genes that belong to another pathway.
#'
#' @param M Binary gene-by-pathway membership matrix.
#' @param n_se Integer. Number of SE genes.
#' @param n_hot Character or integer. Number of hot pathways (see
#'   \code{resolve_n_hot}) (default: 2).
#' @param hot_pathways Character or integer vector, or NULL. Explicit hot
#'   pathways (default: NULL).
#' @param strict Logical. Use only genes that belong exclusively to hot
#'   pathways (default: FALSE).
#' @param prop_outside Numeric in \eqn{[0, 1)}. Proportion of SE genes placed
#'   outside the hot pathways (default: 0).
#'
#' @return A named list with \code{Z} (1 = SE), \code{hot},
#'   \code{n_in}, \code{n_out} and \code{n_shared_used}.
#'
#' @keywords internal
place_se_hot_pathways <- function(M, n_se, n_hot = 2, hot_pathways = NULL, strict = FALSE, prop_outside = 0) {
  M <- as.matrix(M)
  P <- nrow(M)
  K <- ncol(M)
  if (is.null(hot_pathways)) {
    hot <- sort(sample(K, resolve_n_hot(n_hot, K)))
  } else {
    hot <- if (is.character(hot_pathways)) match(hot_pathways, colnames(M)) else as.integer(hot_pathways)
    if (anyNA(hot) || any(hot < 1 | hot > K)) {
      stop("hot_pathways must be pathway names or column numbers of pathway_membership.", call. = FALSE)
    }
    hot <- sort(unique(hot))
  }
  if (prop_outside < 0 || prop_outside >= 1) stop("prop_outside must be in [0, 1).", call. = FALSE)
  n_out <- round(prop_outside * n_se)
  n_in <- n_se - n_out
  in_hot <- rowSums(M[, hot, drop = FALSE]) > 0
  in_cold <- if (length(hot) < K) rowSums(M[, -hot, drop = FALSE]) > 0 else rep(FALSE, P)
  only_hot <- which(in_hot & !in_cold)
  shared <- which(in_hot & in_cold)
  available <- length(only_hot) + if (strict) 0L else length(shared)
  if (available < n_in) {
    stop(sprintf(
      paste0(
        "The %d hot pathway(s) have only %d usable genes (%d only in hot pathways, %d shared ",
        "with other pathways%s), but %d SE genes must go there. Use fewer SE genes (n_se / prop_se), ",
        "more hot pathways (n_hot), fewer and larger pathways (num_pathways)%s."
      ),
      length(hot), available, length(only_hot), length(shared),
      if (strict) ", not allowed with strict = TRUE" else "",
      n_in, if (strict) ", or strict = FALSE" else ""
    ), call. = FALSE)
  }
  hot_order <- if (length(hot) > 1) sample(hot) else hot
  # Take the hot pathways in turn, adding one gene from `pool` to each, until `need` genes are chosen.
  take_in_turn <- function(pool, chosen, need) {
    repeat{
      added <- FALSE
      for (k in hot_order) {
        if (length(chosen) >= need) break
        cand <- setdiff(pool[M[pool, k] == 1], chosen)
        if (length(cand) > 0) {
          chosen <- c(chosen, cand[sample.int(length(cand), 1)])
          added <- TRUE
        }
      }
      if (length(chosen) >= need || !added) {
        return(chosen)
      }
    }
  }
  chosen <- take_in_turn(only_hot, integer(0), n_in)
  if (length(chosen) < n_in && !strict) chosen <- take_in_turn(shared, chosen, n_in)
  if (length(chosen) < n_in) stop("Could not place all SE genes in the hot pathways.", call. = FALSE)
  outside <- integer(0)
  if (n_out > 0) {
    pool_out <- which(!in_hot)
    if (length(pool_out) < n_out) stop("Not enough genes outside the hot pathways for prop_outside.", call. = FALSE)
    outside <- pool_out[sample.int(length(pool_out), n_out)]
  }
  Z <- integer(P)
  Z[c(chosen, outside)] <- 1L
  list(Z = Z, hot = hot, n_in = n_in, n_out = n_out, n_shared_used = sum(chosen %in% shared))
}


## ------------------------------------------------------------
## Function to build a fixed simulation design
## ------------------------------------------------------------

#' Make a Fixed SPHERE Simulation Design
#'
#' Builds a fixed simulation design (spots, library sizes, pathway network,
#' SE status and signal settings) so that every replicate of a scenario is
#' compared against the same truth. Use with \code{\link{simulate_counts}}:
#' call \code{make_sim_design()} once per scenario, then
#' \code{simulate_counts()} once per replicate.
#'
#' @inheritParams generate_genedata_model
#' @param pathway_membership Binary gene-by-pathway membership matrix.
#' @param se_design Character. How SE genes are chosen: \code{"car"}
#'   (pathway-informed through a CAR prior), \code{"hot_pathways"},
#'   \code{"random"}, \code{"exclusive_pathways"} or \code{"exclupath"}
#'   (default: \code{"car"}).
#' @param fix_Z Logical. If TRUE, SE status and signal settings are fixed
#'   across replicates; if FALSE they are redrawn in each replicate. Must be
#'   TRUE for the exclusive designs (default: FALSE).
#' @param n_se Integer or NULL. Number of SE genes for \code{"hot_pathways"}
#'   and \code{"random"} (default: NULL).
#' @param prop_se Numeric. Proportion of SE genes, used when \code{n_se} is
#'   NULL (default: 0.2).
#' @param n_hot Character or integer. Number of hot pathways for
#'   \code{"hot_pathways"} (default: 2).
#' @param hot_pathways Character or integer vector, or NULL. Explicit hot
#'   pathways (default: NULL).
#' @param strict Logical. See \code{place_se_hot_pathways} (default: FALSE).
#' @param prop_outside Numeric in \eqn{[0, 1)}. Proportion of SE genes outside
#'   the hot pathways (default: 0).
#' @param true_Z Integer vector or NULL. True SE indicators for the exclusive
#'   designs (default: NULL).
#' @param se_pathways Integer or character vector, or NULL. SE pathways for
#'   the exclusive designs (default: NULL).
#' @param design_seed Integer or NULL. Seed for the fixed design
#'   (default: NULL).
#'
#' @return A list of class-like type \code{"sphere_design"} containing the
#'   fixed spots, library sizes, SE status, pathway network, gene effects,
#'   signal settings and \code{param}, a list of all design settings.
#'
#' @importFrom stats rnorm median
#' @export
make_sim_design <- function(
    num_spots, num_genes, pathway_membership,
    se_design = c("car", "hot_pathways", "random", "exclusive_pathways", "exclupath"), fix_Z = FALSE,
    base_se_prob = 0.2, rho_u = 0.9, sigma_u = 1.5, car_type = c("leroux", "proper"), standardize_u = TRUE,
    n_se = NULL, prop_se = 0.2, n_hot = 2, hot_pathways = NULL, strict = FALSE, prop_outside = 0,
    true_Z = NULL, se_pathways = NULL,
    sigma_beta = 0.3, eps_sd = 0.25, signal_scenario = "mixed",
    weak_range = c(0.10, 0.25), moderate_range = c(0.25, 0.45), strong_range = c(0.45, 0.70),
    custom_range = c(0.4, 0.7), signal_probs = c(weak = 0.4, moderate = 0.3, strong = 0.3),
    ell_range = c(0.15, 0.35), depth_median = 5000, depth_log_sd = 0.35, baseline_mean = 5,
    x1 = 0, x2 = 1, y1 = 0, y2 = 1, nugget = 1e-8, design_seed = NULL, verbose = TRUE
) {
  se_design <- match.arg(se_design)
  car_type <- match.arg(car_type)

  if (!is.logical(fix_Z) || length(fix_Z) != 1 || is.na(fix_Z)) {
    stop("fix_Z must be TRUE or FALSE.", call. = FALSE)
  }

  # Exclusive pathway membership and SE status must remain fixed together.
  if (se_design %in% c("exclusive_pathways", "exclupath")) {
    if (!isTRUE(fix_Z)) stop("exclusive_pathways requires fix_Z=TRUE.", call. = FALSE)
    if (is.null(true_Z)) stop("Provide true_Z from generate_exclusive_pathway_membership().", call. = FALSE)
    if (is.null(se_pathways)) stop("Provide se_pathways from generate_exclusive_pathway_membership().", call. = FALSE)
  }

  if (!is.null(design_seed)) set.seed(design_seed)

  n <- num_spots
  P <- num_genes
  M <- as.matrix(pathway_membership)

  if (nrow(M) != P) stop("pathway_membership must have num_genes rows.", call. = FALSE)

  # Spots and sequencing depth (same order of random draws as generate_genedata_model).
  spots <- generate_spatial_spots(n, x1, x2, y1, y2)
  domain_scale <- sqrt((x2 - x1) * (y2 - y1))
  ell_range_scaled <- ell_range * domain_scale
  N_i <- round(exp(rnorm(n, log(depth_median), depth_log_sd)))
  N_i <- pmax(N_i, 1L)
  net <- build_pathway_network(P, pathway_membership = M)

  # Which genes are SE.
  alpha0 <- NA_real_
  u <- u_centered <- pi <- rep(NA_real_, P)
  nbr_corr <- NA_real_
  sd_logit_pi <- NA_real_
  hot <- integer(0)
  exclusive_se_pathways <- integer(0)
  n_target <- NA_integer_

  if (se_design == "car") {
    status <- assign_pathway_se_status(P, NULL, M, base_se_prob, rho_u, sigma_u, car_type = car_type, standardize = standardize_u)
    Z <- status$Z
    alpha0 <- status$alpha0
    u <- status$u
    u_centered <- status$u_centered
    pi <- status$pi
    nbr_corr <- status$nbr_corr
    sd_logit_pi <- status$sd_logit_pi
  } else if (se_design == "exclupath") {
    Z <- as.integer(true_Z)
    if (length(Z) != P || anyNA(Z) || any(!Z %in% c(0, 1))) stop("true_Z must be a binary vector of length num_genes.", call. = FALSE)
    check <- validate_exclupath(M, Z, se_pathways)
    if (!isTRUE(check$valid)) stop(paste("Invalid exclupath structure:", paste(names(check$checks)[!check$checks], collapse = ", ")), call. = FALSE)
    exclusive_se_pathways <- sort(unique(as.integer(se_pathways)))
    n_target <- sum(Z)
    exclupath_se_per <- median(colSums(M[Z == 1, exclusive_se_pathways, drop = FALSE]))
    exclupath_nonse_per <- median(colSums(M[Z == 0, exclusive_se_pathways, drop = FALSE]))
  } else if (se_design == "exclusive_pathways") {
    # Use the SE indicators generated jointly with the pathway membership.
    Z <- as.integer(true_Z)

    if (length(Z) != P || anyNA(Z) || any(!Z %in% c(0, 1))) {
      stop("true_Z must be a binary vector of length num_genes.", call. = FALSE)
    }

    # Verify that SE and non-SE genes belong exclusively to their designated pathways.
    check <- validate_exclusive_pathways(M, Z, se_pathways)

    if (!isTRUE(check$valid)) {
      stop(paste("Invalid exclusive pathway structure:", paste(names(check$checks)[!check$checks], collapse = ", ")), call. = FALSE)
    }

    exclusive_se_pathways <- if (is.character(se_pathways)) {
      match(se_pathways, colnames(M))
    } else {
      as.integer(se_pathways)
    }

    exclusive_se_pathways <- sort(unique(exclusive_se_pathways))
    n_target <- sum(Z)

    if (!is.null(n_se) && n_se != n_target) {
      stop("n_se does not match the number of SE genes in true_Z.", call. = FALSE)
    }
  } else {
    n_target <- as.integer(if (is.null(n_se)) round(prop_se * P) else n_se)

    if (length(n_target) != 1 || is.na(n_target) || n_target < 1 || n_target >= P) {
      stop("n_se must be between 1 and num_genes-1.", call. = FALSE)
    }

    if (se_design == "hot_pathways") {
      pl <- place_se_hot_pathways(M, n_target, n_hot, hot_pathways, strict, prop_outside)
      Z <- pl$Z
      hot <- pl$hot
    } else {
      Z <- integer(P)
      Z[sample.int(P, n_target)] <- 1L
    }
  }

  se_idx <- which(Z == 1)
  nonse_idx <- which(Z == 0)

  # Gene levels and SE signal strength / width (fixed when fix_Z=TRUE).
  Beta <- rnorm(P, 0, sigma_beta)
  gp <- assign_gp_parameters(
    P = P, se_idx = se_idx, signal_scenario = signal_scenario, sigma_eta_range = custom_range,
    weak_range = weak_range, moderate_range = moderate_range, strong_range = strong_range,
    signal_probs = signal_probs, ell_range = ell_range_scaled
  )

  # Gene names remain the same even when SE status changes.
  gene_names <- paste0("Gene", seq_len(P))

  design <- list(
    type = "sphere_design",
    fix_Z = fix_Z,
    spots = spots, N_i = N_i, gene_names = gene_names,
    Z = Z, se_idx = se_idx, nonse_idx = nonse_idx,
    alpha0 = alpha0, u = u, u_centered = u_centered, pi = pi, nbr_corr = nbr_corr, sd_logit_pi = sd_logit_pi,
    hot_pathways = if (length(hot)) colnames(M)[hot] else character(0),
    hot_pathway_idx = hot,
    se_pathways = exclusive_se_pathways,
    W = net$W, degree = net$degree, pathway_membership = M, shared_pathways = net$shared_pathways,
    Beta = Beta, sigma_eta = gp$sigma_eta, ell = gp$ell, signal_class = gp$signal_class,
    domain_scale = domain_scale, ell_range_scaled = ell_range_scaled,
    param = list(
      n = n, p = P, se_design = se_design, fix_Z = fix_Z,
      design_seed = if (is.null(design_seed)) NA else design_seed,
      base_se_prob = if (se_design == "car") base_se_prob else NA,
      rho_u = if (se_design == "car") rho_u else NA,
      sigma_u = if (se_design == "car") sigma_u else NA,
      car_type = if (se_design == "car") car_type else NA,
      standardize_u = standardize_u,
      n_se = sum(Z), n_se_target = n_target,
      n_hot = if (se_design == "hot_pathways") length(hot) else NA,
      n_se_pathways = if (se_design %in% c("exclusive_pathways", "exclupath")) length(exclusive_se_pathways) else NA,
      se_per_pathway = if (se_design == "exclupath") exclupath_se_per else NA,
      nonse_per_se_pathway = if (se_design == "exclupath") exclupath_nonse_per else NA,
      strict = strict, prop_outside = if (se_design == "hot_pathways") prop_outside else NA,
      n_pathways = ncol(M), sigma_beta = sigma_beta, eps_sd = eps_sd,
      signal_scenario = signal_scenario, weak_range = weak_range, moderate_range = moderate_range,
      strong_range = strong_range, custom_range = custom_range,
      signal_probs = signal_probs, ell_range = ell_range,
      domain_scale = domain_scale, ell_range_scaled = ell_range_scaled, depth_median = depth_median,
      baseline_mean = baseline_mean, nugget = nugget
    )
  )

  if (verbose) {
    cat("---- SPHERE simulation design:", design_label(design$param), "----\n")
    cat("Fix Z:", fix_Z, "\n")
    cat("Spots:", n, " Genes:", P, " Pathways:", ncol(M), " Mean network degree:", round(mean(net$degree), 2), "\n")
    if (length(hot)) cat("Hot pathways:", paste(design$hot_pathways, collapse = ", "), "\n")
    if (se_design %in% c("exclusive_pathways", "exclupath")) {
      cat("SE-only pathways:", paste(colnames(M)[exclusive_se_pathways], collapse = ", "), "\n")
      cat("Exclusive structure valid:", TRUE, "\n")
    }
    if (sum(Z) > 0) cat("Signal classes:", paste(names(table(gp$signal_class[se_idx])), table(gp$signal_class[se_idx]), collapse = ", "), "\n")
    invisible(summarize_pathway_se(design))
  }
  design
}


## ------------------------------------------------------------
## Function to draw one spatial pattern shape
## ------------------------------------------------------------

#' Draw a Random Spatial Pattern Shape
#'
#' Draws a random location and orientation for one pattern, in coordinates
#' rescaled to \eqn{[0, 1]}. Centres stay away from the edges so hotspots and
#' rings are not cut off.
#'
#' @return A named list with \code{cx}, \code{cy}, \code{angle},
#'   \code{radius}, \code{freq}, \code{phase} and \code{wave_amp}.
#'
#' @importFrom stats runif
#' @keywords internal
draw_pattern_shape <- function() {
  list(cx = runif(1, 0.25, 0.75), cy = runif(1, 0.25, 0.75),   # centre
       angle = runif(1, 0, 2 * pi),                              # direction (streak, gradient, wave); full circle so gradients can point any way
       radius = runif(1, 0.20, 0.30),                            # ring radius
       freq = runif(1, 1, 2), phase = runif(1, 0, 2 * pi),      # wave: cycles across the slide, offset
       wave_amp = runif(1, 0.10, 0.20))                          # wave: height of the curve
}


## ------------------------------------------------------------
## Function to simulate a named spatial pattern
## ------------------------------------------------------------

#' Simulate a Spatial Expression Pattern
#'
#' Returns a spatial signal \eqn{\eta \in [0, 1]} for every spot
#' (1 = high-expression region, 0 = background):
#' \describe{
#'   \item{hotspot}{Disc of the \code{high_prop} spots closest to the centre.}
#'   \item{streak}{Band of the \code{high_prop} spots closest to a line
#'     through the centre (any angle).}
#'   \item{ring}{The \code{high_prop} spots closest to a circle around the
#'     centre.}
#'   \item{wave}{The \code{high_prop} spots closest to a sine curve (any
#'     orientation).}
#'   \item{gradient}{Linear ramp along a random direction: 0 before
#'     \code{grad_start}, rising to 1.}
#' }
#'
#' @param spots A data frame with columns \code{x} and \code{y}.
#' @param pat_type Character. One of \code{"hotspot"}, \code{"streak"},
#'   \code{"gradient"}, \code{"ring"} or \code{"wave"}.
#' @param high_prop Numeric. Proportion of spots in the high-expression
#'   region (default: 0.3).
#' @param grad_start Numeric. Start of the gradient ramp in rescaled
#'   coordinates (default: 0.3).
#' @param shape List or NULL. A shape from \code{draw_pattern_shape}; reuse it
#'   so several genes share one pattern. If NULL, a new shape is drawn
#'   (default: NULL).
#'
#' @return A named list with \code{eta}, \code{mark} (\code{"high"},
#'   \code{"medium"} or \code{"low"}) and the \code{shape} used.
#'
#' @export
simulate_pattern <- function(spots, pat_type = c("hotspot", "streak", "gradient", "ring", "wave"),
                             high_prop = 0.3, grad_start = 0.3, shape = NULL) {
  pat_type <- match.arg(pat_type)
  xy <- as.matrix(spots[, c("x", "y")])
  u <- (xy[, 1] - min(xy[, 1])) / diff(range(xy[, 1]))
  v <- (xy[, 2] - min(xy[, 2])) / diff(range(xy[, 2]))
  n <- nrow(xy)
  n_high <- max(1L, round(high_prop * n))
  if (is.null(shape)) shape <- draw_pattern_shape()
  along  <-  (u - shape$cx) * cos(shape$angle) + (v - shape$cy) * sin(shape$angle)
  across <- -(u - shape$cx) * sin(shape$angle) + (v - shape$cy) * cos(shape$angle)
  closest <- function(score) { e <- numeric(n); e[order(score)[seq_len(n_high)]] <- 1; e }
  eta <- switch(pat_type,
                hotspot  = closest((u - shape$cx)^2 + (v - shape$cy)^2),
                streak   = closest(abs(across)),
                ring     = closest(abs(sqrt((u - shape$cx)^2 + (v - shape$cy)^2) - shape$radius)),
                wave     = closest(abs(across - shape$wave_amp * sin(2 * pi * shape$freq * along + shape$phase))),
                gradient = { t <- (along - min(along)) / diff(range(along)); pmin(1, pmax(0, (t - grad_start) / (1 - grad_start))) })
  mark <- if (pat_type == "gradient") ifelse(eta >= 0.75, "high", ifelse(eta <= 0.40, "low", "medium")) else
    ifelse(eta >= 0.5, "high", "low")
  list(eta = eta, mark = mark, shape = shape)
}


## ------------------------------------------------------------
## Function to assign pattern type, strength and size for SE genes
## ------------------------------------------------------------

#' Assign Pattern Type, Strength and Size for SE Genes
#'
#' Assigns each SE gene a pattern shape, a strength (log fold change between
#' the high region and the background) and a high-region size. Strength
#' classes weak / moderate / strong are drawn exactly like the model-based
#' signal classes.
#'
#' @param P Integer. Number of genes.
#' @param se_idx Integer vector. Indices of SE genes.
#' @param pattern_types Character. \code{"mixed"} (all five shapes in turn)
#'   or a subset such as \code{c("hotspot", "streak")} (default:
#'   \code{"mixed"}).
#' @param signal_scenario Character. One of \code{"weak"},
#'   \code{"moderate"}, \code{"strong"}, \code{"mixed"} or \code{"custom"}
#'   (default: \code{"mixed"}).
#' @param weak_range Numeric vector of length 2. Log fold change range for
#'   weak signals (default: \code{c(0.5, 0.8)}).
#' @param moderate_range Numeric vector of length 2. Log fold change range
#'   for moderate signals (default: \code{c(0.8, 1.1)}).
#' @param strong_range Numeric vector of length 2. Log fold change range for
#'   strong signals (default: \code{c(1.1, 1.5)}).
#' @param custom_range Numeric vector of length 2. Log fold change range for
#'   \code{"custom"} (default: \code{c(0.8, 1.2)}).
#' @param signal_probs Named numeric vector (weak, moderate, strong) summing
#'   to 1 (default: \code{c(weak = 0.4, moderate = 0.3, strong = 0.3)}).
#' @param high_prop_range Numeric vector of length 2. Range of the
#'   high-region proportion (default: \code{c(0.2, 0.4)}).
#' @param pathway_membership Binary matrix or NULL. Needed when
#'   \code{pathway_patterns = "shared"} (default: NULL).
#' @param pathway_patterns Character. \code{"independent"}: every SE gene has
#'   its own shape and location; \code{"shared"}: SE genes of the same SE
#'   pathway share one shape, location and high-region size (default:
#'   \code{"independent"}).
#'
#' @return A named list with \code{pattern_type}, \code{pattern_amp},
#'   \code{pattern_prop}, \code{pattern_module} and \code{signal_class}.
#'
#' @importFrom stats runif setNames
#' @keywords internal
assign_pattern_parameters <- function(P, se_idx, pattern_types = "mixed", signal_scenario = "mixed",
                                      weak_range = c(0.5, 0.8), moderate_range = c(0.8, 1.1),
                                      strong_range = c(1.1, 1.5), custom_range = c(0.8, 1.2),
                                      signal_probs = c(weak = 0.4, moderate = 0.3, strong = 0.3),
                                      high_prop_range = c(0.2, 0.4), pathway_membership = NULL,
                                      pathway_patterns = c("independent", "shared")) {
  pathway_patterns <- match.arg(pathway_patterns)
  available <- c("hotspot", "streak", "gradient", "ring", "wave")
  types <- if (identical(pattern_types, "mixed")) available else match.arg(pattern_types, available, several.ok = TRUE)
  if (pathway_patterns == "shared" && is.null(pathway_membership))
    stop("pathway_patterns = 'shared' needs pathway_membership.", call. = FALSE)
  # Strength classes and values: same class sampling as the model-based signal.
  st <- assign_gp_parameters(P = P, se_idx = se_idx, signal_scenario = signal_scenario,
                             sigma_eta_range = custom_range, weak_range = weak_range,
                             moderate_range = moderate_range, strong_range = strong_range,
                             signal_probs = signal_probs, ell_range = c(1, 1))
  amp <- st$sigma_eta
  module <- rep(NA_integer_, P)
  if (length(se_idx)) {
    if (pathway_patterns == "shared") {
      M <- as.matrix(pathway_membership)
      n_se_pw <- colSums(M[se_idx, , drop = FALSE])
      module[se_idx] <- vapply(se_idx, function(j) {
        pw <- which(M[j, ] == 1)
        if (!length(pw)) return(-j)                         # gene without a pathway: its own module
        pw[which.max(n_se_pw[pw])]
      }, integer(1))
    } else module[se_idx] <- se_idx
  }
  mods <- sort(unique(module[se_idx]))
  # One shape type and one high-region size per module, shapes dealt out in turn so all are used.
  mod_type <- setNames(if (length(mods)) rep_len(types[sample.int(length(types))], length(mods)) else character(0), mods)
  mod_prop <- setNames(runif(length(mods), high_prop_range[1], high_prop_range[2]), mods)
  type <- rep("NonSE", P); prop <- rep(NA_real_, P)
  type[se_idx] <- mod_type[as.character(module[se_idx])]
  prop[se_idx] <- mod_prop[as.character(module[se_idx])]
  list(pattern_type = type, pattern_amp = amp, pattern_prop = prop, pattern_module = module,
       signal_class = st$signal_class)
}


## ------------------------------------------------------------
## Function to build a fixed pattern-based design
## ------------------------------------------------------------

#' Make a Fixed Pattern-Based SPHERE Design
#'
#' Builds the same design as \code{\link{make_sim_design}} (which genes are
#' SE and which pathways they sit in), but replaces the GP signal of each SE
#' gene with a fixed spatial shape (hotspot, streak, gradient, ring, wave).
#' Counts are generated the same way. Use with \code{\link{simulate_counts}}.
#'
#' @param ... Arguments passed to \code{\link{make_sim_design}}.
#' @param pattern_types Character. \code{"mixed"} or a subset of the five
#'   shapes (default: \code{"mixed"}).
#' @param pathway_patterns Character. \code{"independent"} or \code{"shared"}
#'   (default: \code{"independent"}).
#' @param pattern_signal_scenario Character or NULL. Strength scenario for the
#'   patterns. If NULL, uses the design's \code{signal_scenario}
#'   (default: NULL).
#' @param pattern_weak_range Numeric vector of length 2. Log fold change range
#'   for weak patterns (default: \code{c(0.5, 0.8)}).
#' @param pattern_moderate_range Numeric vector of length 2. Log fold change
#'   range for moderate patterns (default: \code{c(0.8, 1.1)}).
#' @param pattern_strong_range Numeric vector of length 2. Log fold change
#'   range for strong patterns (default: \code{c(1.1, 1.5)}).
#' @param pattern_custom_range Numeric vector of length 2. Log fold change
#'   range for \code{"custom"} (default: \code{c(0.8, 1.2)}).
#' @param pattern_signal_probs Named numeric vector or NULL. Class
#'   probabilities. If NULL, uses the design's \code{signal_probs}
#'   (default: NULL).
#' @param high_prop_range Numeric vector of length 2. Range of the
#'   high-region proportion (default: \code{c(0.2, 0.4)}).
#' @param grad_start Numeric. Start of the gradient ramp (default: 0.3).
#' @param center_pattern Logical. If TRUE, \eqn{\eta} is centred so a pattern
#'   moves expression around the slide without raising the gene's overall
#'   mean (default: TRUE).
#' @param fix_pattern_shape Logical. If TRUE, the same shapes are used in
#'   every replicate; if FALSE, a new location / orientation is drawn each
#'   replicate (default: FALSE).
#' @param zero_inflation Numeric. Extra share of counts set to 0
#'   (default: 0).
#'
#' @return A \code{"sphere_design"} list as from \code{make_sim_design}, with
#'   pattern settings added and \code{param$signal_model = "pattern"}.
#'
#' @importFrom stats setNames na.omit
#' @export
make_pattern_design <- function(..., pattern_types = "mixed", pathway_patterns = c("independent", "shared"),
                                pattern_signal_scenario = NULL,
                                pattern_weak_range = c(0.5, 0.8), pattern_moderate_range = c(0.8, 1.1),
                                pattern_strong_range = c(1.1, 1.5), pattern_custom_range = c(0.8, 1.2),
                                pattern_signal_probs = NULL, high_prop_range = c(0.2, 0.4), grad_start = 0.3,
                                center_pattern = TRUE, fix_pattern_shape = FALSE, zero_inflation = 0) {
  pathway_patterns <- match.arg(pathway_patterns)
  args <- list(...)
  verbose <- if (is.null(args$verbose)) TRUE else isTRUE(args$verbose)
  args$verbose <- FALSE
  design <- do.call(make_sim_design, args)
  pr <- design$param
  P <- length(design$Z)
  scenario <- if (is.null(pattern_signal_scenario)) pr$signal_scenario else pattern_signal_scenario
  probs <- if (is.null(pattern_signal_probs)) pr$signal_probs else pattern_signal_probs
  pp <- assign_pattern_parameters(P, design$se_idx, pattern_types, scenario,
                                  pattern_weak_range, pattern_moderate_range, pattern_strong_range,
                                  pattern_custom_range, probs, high_prop_range,
                                  design$pathway_membership, pathway_patterns)
  design$signal_class <- pp$signal_class
  design$sigma_eta <- pp$pattern_amp          # pattern strength (kept under this name for summaries)
  design$ell <- rep(NA_real_, P)
  design[c("pattern_type", "pattern_amp", "pattern_prop", "pattern_module")] <- pp[c("pattern_type", "pattern_amp", "pattern_prop", "pattern_module")]
  mods <- sort(unique(stats::na.omit(pp$pattern_module)))
  design$pattern_shapes <- if (fix_pattern_shape) setNames(lapply(mods, function(m) draw_pattern_shape()), mods) else NULL
  drop <- c("weak_range", "moderate_range", "strong_range", "custom_range", "ell_range", "ell_range_scaled", "nugget",
            "signal_scenario", "signal_probs")   # replaced below; dropping them avoids two copies of the same setting
  design$param <- c(pr[setdiff(names(pr), drop)], list(
    signal_model = "pattern", signal_scenario = scenario, signal_probs = probs,
    pattern_types = pattern_types, pathway_patterns = pathway_patterns,
    weak_range = pattern_weak_range, moderate_range = pattern_moderate_range,
    strong_range = pattern_strong_range, custom_range = pattern_custom_range,
    high_prop_range = high_prop_range, grad_start = grad_start, center_pattern = center_pattern,
    fix_pattern_shape = fix_pattern_shape, zero_inflation = zero_inflation))
  if (verbose) {
    M <- design$pathway_membership
    cat("---- Pattern-based SPHERE design:", design_label(design$param), "----\n")
    cat("Fix Z:", design$fix_Z, " | Spots:", nrow(design$spots), " Genes:", P, " Pathways:", ncol(M), "\n")
    if (length(design$se_pathways)) cat("SE pathways:", paste(colnames(M)[design$se_pathways], collapse = ", "), "\n")
    se <- design$se_idx
    if (length(se)) {
      cat("Pattern types:", paste(names(table(pp$pattern_type[se])), table(pp$pattern_type[se]), collapse = ", "), "\n")
      cat("Signal classes:", paste(names(table(pp$signal_class[se])), table(pp$signal_class[se]), collapse = ", "), "\n")
      cat("Log fold change (high vs background):", paste(round(range(pp$pattern_amp[se]), 2), collapse = " - "), "\n")
      cat("Pathway patterns:", pathway_patterns, "(", length(mods), "pattern modules )\n")
    }
    invisible(summarize_pathway_se(design))
  }
  design
}


## ------------------------------------------------------------
## Function to simulate one replicate from a design
## ------------------------------------------------------------

#' Simulate One Replicate from a SPHERE Design
#'
#' Generates one count dataset from a design made by
#' \code{\link{make_sim_design}} or \code{\link{make_pattern_design}}. When
#' \code{fix_Z = FALSE}, SE status, gene effects and signal settings are
#' redrawn; otherwise only the spatial patterns, noise and counts are redrawn.
#' Exclusive pathway designs always keep their fixed membership and SE status.
#'
#' @param design A \code{"sphere_design"} list.
#' @param seed Integer or NULL. Random seed for this replicate (default: NULL).
#' @param rep_id Integer. Replicate index stored in \code{param}
#'   (default: NA).
#' @param verbose Logical. If TRUE, prints a one-line summary
#'   (default: FALSE).
#'
#' @return A named list with \code{Y}, \code{mu0}, \code{epsilon}, \code{eta},
#'   \code{loglambda}, \code{lambda}, \code{Y_mean}, the design quantities for
#'   this replicate and \code{param} (design settings plus \code{seed} and
#'   \code{rep_id}).
#'
#' @importFrom stats rnorm rpois runif setNames na.omit
#' @export
simulate_counts <- function(design, seed = NULL, rep_id = NA, verbose = FALSE) {
  if (!identical(design$type, "sphere_design")) stop("design must come from make_sim_design() or make_pattern_design().", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  pr <- design$param
  is_pattern <- identical(pr$signal_model, "pattern")
  n <- nrow(design$spots)
  P <- length(design$Z)
  coords <- as.matrix(design$spots[, c("x", "y")])

  if (pr$se_design %in% c("exclusive_pathways", "exclupath") && !isTRUE(design$fix_Z)) {
    stop("exclusive_pathways requires fix_Z=TRUE.", call. = FALSE)
  }

  d <- design

  if (!isTRUE(design$fix_Z)) {
    if (pr$se_design == "car") {
      status <- assign_pathway_se_status(
        P, NULL, design$pathway_membership,
        pr$base_se_prob, pr$rho_u, pr$sigma_u,
        car_type = pr$car_type, standardize = pr$standardize_u
      )
      d$Z <- status$Z
      d$alpha0 <- status$alpha0
      d$u <- status$u
      d$u_centered <- status$u_centered
      d$pi <- status$pi
      d$nbr_corr <- status$nbr_corr
      d$sd_logit_pi <- status$sd_logit_pi
    } else if (pr$se_design == "hot_pathways") {
      pl <- place_se_hot_pathways(
        design$pathway_membership, pr$n_se_target,
        n_hot = pr$n_hot,
        hot_pathways = design$hot_pathway_idx,
        strict = pr$strict, prop_outside = pr$prop_outside
      )
      d$Z <- pl$Z
    } else {
      d$Z <- integer(P)
      d$Z[sample.int(P, pr$n_se_target)] <- 1L
    }

    d$se_idx <- which(d$Z == 1)
    d$nonse_idx <- which(d$Z == 0)

    d$Beta <- rnorm(P, 0, pr$sigma_beta)
    if (is_pattern) {
      pp <- assign_pattern_parameters(P, d$se_idx, pr$pattern_types, pr$signal_scenario,
                                      pr$weak_range, pr$moderate_range, pr$strong_range, pr$custom_range,
                                      pr$signal_probs, pr$high_prop_range, design$pathway_membership,
                                      pr$pathway_patterns)
      d[c("pattern_type", "pattern_amp", "pattern_prop", "pattern_module")] <- pp[c("pattern_type", "pattern_amp", "pattern_prop", "pattern_module")]
      d$sigma_eta <- pp$pattern_amp
      d$signal_class <- pp$signal_class
      d$pattern_shapes <- NULL                 # new SE genes -> new shapes
    } else {
      gp <- assign_gp_parameters(
        P = P, se_idx = d$se_idx, signal_scenario = pr$signal_scenario,
        sigma_eta_range = pr$custom_range,
        weak_range = pr$weak_range, moderate_range = pr$moderate_range,
        strong_range = pr$strong_range, signal_probs = pr$signal_probs,
        ell_range = pr$ell_range_scaled
      )
      d$sigma_eta <- gp$sigma_eta
      d$ell <- gp$ell
      d$signal_class <- gp$signal_class
    }
  }

  d$param$n_se <- sum(d$Z)

  epsilon <- matrix(rnorm(n * P, 0, pr$eps_sd), n, P)
  eta <- matrix(0, n, P)

  if (is_pattern) {
    # One shape per pattern module (one gene, or one shared pathway module); a new location in each
    # replicate unless the design fixed the shapes.
    mods <- sort(unique(stats::na.omit(d$pattern_module[d$se_idx])))
    shapes <- if (!is.null(d$pattern_shapes)) d$pattern_shapes else setNames(lapply(mods, function(m) draw_pattern_shape()), mods)
    for (j in d$se_idx) {
      m <- as.character(d$pattern_module[j])
      e <- simulate_pattern(design$spots, d$pattern_type[j], d$pattern_prop[j], pr$grad_start, shapes[[m]])$eta
      eta[, j] <- d$pattern_amp[j] * (e - if (isTRUE(pr$center_pattern)) mean(e) else 0)
    }
    d$pattern_shapes <- shapes
  } else {
    for (j in d$se_idx) eta[, j] <- simulate_full_gp(coords, d$sigma_eta[j], d$ell[j], pr$nugget)
  }

  mu0 <- log(pr$baseline_mean / pr$depth_median)
  loglambda <- matrix(mu0, n, P) + matrix(d$Beta, n, P, byrow = TRUE) + epsilon + eta
  lambda <- exp(loglambda)
  Y_mean <- lambda * d$N_i
  Y <- matrix(rpois(n * P, as.vector(Y_mean)), n, P)
  if (is_pattern && isTRUE(pr$zero_inflation > 0)) Y[matrix(runif(n * P) < pr$zero_inflation, n, P)] <- 0L
  gene_names <- character(P)
  gene_names[d$nonse_idx] <- paste0("Gene", d$nonse_idx, "_NonSE")
  gene_names[d$se_idx] <- paste0("Gene", d$se_idx, "_SE_", d$signal_class[d$se_idx])
  colnames(Y) <- gene_names
  colnames(loglambda) <- gene_names
  colnames(eta) <- gene_names

  if (verbose) cat("Replicate", rep_id, "| seed", seed, "| SE genes", sum(d$Z), "| mean count", round(mean(Y), 3), "\n")

  keep <- c(
    "spots", "N_i", "Beta", "Z", "se_idx", "nonse_idx", "alpha0", "u", "u_centered", "pi", "W", "degree",
    "pathway_membership", "shared_pathways", "sigma_eta", "ell", "signal_class", "domain_scale",
    "ell_range_scaled", "nbr_corr", "sd_logit_pi", "hot_pathways", "gene_names", "fix_Z", "se_pathways",
    if (is_pattern) c("pattern_type", "pattern_amp", "pattern_prop", "pattern_module", "pattern_shapes")
  )

  c(
    list(Y = Y, mu0 = mu0, epsilon = epsilon, eta = eta, loglambda = loglambda, lambda = lambda, Y_mean = Y_mean),
    d[intersect(keep, names(d))],
    list(param = c(d$param, list(seed = if (is.null(seed)) NA else seed, rep_id = rep_id)))
  )
}
