## ============================================================
## SPHERE: Helper Functions
## ============================================================


## ------------------------------------------------------------
## Color palette for spatial expression plots
## ------------------------------------------------------------

pal <- colorRampPalette(c('#00274c', '#00274c', "lightyellow2",
                          '#ffcb05', '#ffcb05'))

#' Plot Gene Expression Data
#'
#' Creates spatial expression plots for selected genes arranged in a grid.
#'
#' @param gene_df A matrix or data frame of gene expression values
#'   (spots x genes).
#' @param spot A matrix of spatial coordinates with two columns (x, y).
#' @param genes Character vector of gene names to plot. If NULL, plots
#'   all genes (default: NULL).
#' @param nrow Integer. Number of rows in the plot grid (default: 1).
#'
#' @return A grid of ggplot objects arranged by \pkg{gridExtra}.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#'   theme_minimal labs theme element_blank
#' @export
plot_gene_data <- function(gene_df, spot, genes = NULL, nrow = 1) {
  x <- y <- NULL  # fix R CMD CHECK note
  gene_df <- as.matrix(gene_df)
  if (is.null(genes)) genes <- colnames(gene_df)

  rel_expr <- apply(gene_df, 2, function(x) {
    if (max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  })
  rel_expr <- as.matrix(rel_expr)

  plots <- lapply(genes, function(gene_name) {
    df_plot <- data.frame(x = spot[, 1], y = spot[, 2],
                          expression = rel_expr[, gene_name])
    ggplot(df_plot, aes(x, y, color = expression)) +
      geom_point(size = 4) +
      scale_color_gradientn(colours = pal(5), limits = c(0, 1)) +
      theme_minimal() +
      labs(title = paste("Relative Expression for\n", gene_name),
           color = "Relative\nExpression") +
      theme(axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank())
  })
  do.call(gridExtra::grid.arrange, c(plots, nrow = nrow))
}


## ------------------------------------------------------------
## Function to plot simulated spatial patterns
## ------------------------------------------------------------

#' Plot Simulated Spatial Patterns
#'
#' Spatial map of the signal (\eqn{\eta}) or counts of selected genes in one
#' simulated replicate, one panel per gene, titled with the gene, its pattern,
#' strength class and pathway(s).
#'
#' @param d A simulated dataset from \code{\link{simulate_counts}}.
#' @param genes Integer vector. Genes to plot (default: the SE genes).
#' @param what Character. \code{"eta"} (true signal) or \code{"count"}
#'   (\code{log(1 + count / depth * median depth)}) (default: \code{"eta"}).
#' @param ncol Integer. Number of panel columns (default: 5).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom stats median
#' @export
plot_sim_patterns <- function(d, genes = d$se_idx, what = c("eta", "count"), ncol = 5) {
  x <- y <- value <- gene <- NULL  # fix R CMD CHECK note
  what <- match.arg(what)
  M <- d$pathway_membership
  lab <- vapply(genes, function(j) sprintf("Gene%d | %s, %s | %s", j,
                                           if (!is.null(d$pattern_type)) d$pattern_type[j] else "GP", d$signal_class[j],
                                           paste(colnames(M)[M[j, ] == 1], collapse = "+")), "")
  val <- if (what == "eta") d$eta[, genes, drop = FALSE] else log1p(d$Y[, genes, drop = FALSE] / d$N_i * stats::median(d$N_i))
  df <- data.frame(x = rep(d$spots$x, length(genes)), y = rep(d$spots$y, length(genes)),
                   value = as.vector(val), gene = factor(rep(lab, each = nrow(d$spots)), levels = lab))
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = value)) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_viridis_c(name = if (what == "eta") "signal" else "log count") +
    ggplot2::facet_wrap(~ gene, ncol = ncol) + ggplot2::coord_equal() + ggplot2::theme_bw(base_size = 8) +
    ggplot2::labs(x = NULL, y = NULL)
}


## ------------------------------------------------------------
## Function to build the pathway network
## ------------------------------------------------------------

#' Build the Pathway Network
#'
#' Builds the gene network in which genes sharing at least one pathway are
#' neighbours. Supply either \code{gene_grp} (one pathway per gene) or
#' \code{pathway_membership} (genes may belong to several pathways).
#'
#' @param num_genes Integer. Number of genes.
#' @param gene_grp Vector or NULL. Single pathway label per gene
#'   (default: NULL).
#' @param pathway_membership Binary matrix or NULL. Gene-by-pathway
#'   membership (default: NULL).
#'
#' @return A named list with \code{W} (binary adjacency matrix),
#'   \code{degree}, \code{membership} and \code{shared_pathways}.
#'
#' @keywords internal
build_pathway_network <- function(num_genes,gene_grp=NULL,pathway_membership=NULL){
  if(is.null(gene_grp)&&is.null(pathway_membership)) stop("Provide either gene_grp or pathway_membership.")
  if(!is.null(gene_grp)&&!is.null(pathway_membership)) stop("Provide gene_grp or pathway_membership, not both.")
  if(!is.null(gene_grp)){
    if(length(gene_grp)!=num_genes) stop("gene_grp must have length num_genes.")
    groups <- unique(gene_grp)
    M <- matrix(0,num_genes,length(groups))
    for(k in seq_along(groups)) M[,k] <- as.integer(gene_grp==groups[k])
    colnames(M) <- paste0("Pathway_",groups)
  }else{
    M <- as.matrix(pathway_membership)
    if(nrow(M)!=num_genes) stop("pathway_membership must have num_genes rows.")
    if(any(!M%in%c(0,1))) stop("pathway_membership must contain only 0 and 1.")
    if(any(rowSums(M)==0)) stop("Every gene must belong to at least one pathway.")
  }
  shared <- tcrossprod(M)
  W <- 1L*(shared>0)
  diag(W) <- 0L
  degree <- rowSums(W)
  list(W=W,degree=degree,membership=M,shared_pathways=shared)
}


## ------------------------------------------------------------
## Function to permute the network (wrong biology)
## ------------------------------------------------------------

#' Permute a Pathway Network
#'
#' Shuffles gene labels in an adjacency matrix, keeping its structure but
#' breaking its link to the true biology.
#'
#' @param W Binary symmetric adjacency matrix.
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return The permuted adjacency matrix.
#'
#' @keywords internal
permute_network<-function(W,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  P<-nrow(W)
  perm<-sample(seq_len(P))
  W_perm<-W[perm,perm,drop=FALSE]
  dimnames(W_perm)<-dimnames(W)
  return(W_perm)
}


## ------------------------------------------------------------
## Function to partially misspecify the network
## ------------------------------------------------------------

#' Rewire a Proportion of Network Edges
#'
#' Partially misspecifies a network by moving a proportion of its edges to
#' random gene pairs, keeping the number of edges fixed.
#'
#' @param W Binary symmetric adjacency matrix.
#' @param proportion Numeric in \eqn{[0, 1]}. Proportion of edges to move
#'   (default: 0.2).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return The rewired adjacency matrix.
#'
#' @keywords internal
rewire_network<-function(W,proportion=0.2,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  W<-as.matrix(W)
  if(nrow(W)!=ncol(W))stop("W must be square.")
  if(!isTRUE(all(W==t(W))))stop("W must be symmetric.")
  if(!all(W%in%c(0,1)))stop("W must contain only 0 and 1.")
  if(proportion<0||proportion>1)stop("proportion must be between 0 and 1.")
  diag(W)<-0
  edges<-which(upper.tri(W)&W==1,arr.ind=TRUE)
  n_edges<-nrow(edges)
  n_change<-round(proportion*n_edges)
  if(n_change==0)return(W)
  remove_idx<-sample(seq_len(n_edges),n_change)
  removed<-edges[remove_idx,,drop=FALSE]
  W[cbind(removed[,1],removed[,2])]<-0
  W[cbind(removed[,2],removed[,1])]<-0
  available<-which(upper.tri(W)&W==0,arr.ind=TRUE)
  removed_keys<-paste(removed[,1],removed[,2],sep="_")
  available_keys<-paste(available[,1],available[,2],sep="_")
  available<-available[!(available_keys%in%removed_keys),,drop=FALSE]
  if(nrow(available)<n_change)stop("Not enough new connections available.")
  add_idx<-sample(seq_len(nrow(available)),n_change)
  added<-available[add_idx,,drop=FALSE]
  W[cbind(added[,1],added[,2])]<-1
  W[cbind(added[,2],added[,1])]<-1
  return(W)
}


## ------------------------------------------------------------
## Function to prepare the sparse CAR inputs for Stan
## ------------------------------------------------------------

#' Prepare Sparse CAR Inputs for Stan
#'
#' Converts a binary adjacency matrix into the sparse inputs required by the
#' SPHERE Stan model: connected and isolated gene indices, an edge list,
#' neighbour counts and the eigenvalues used for the Leroux and sparse proper
#' CAR log determinants.
#'
#' @param W Binary symmetric adjacency matrix.
#'
#' @return A named list with \code{P_car}, \code{car_index}, \code{P_iso},
#'   \code{iso_index}, \code{W_n}, \code{W_sparse}, \code{D_sparse},
#'   \code{lap_eig} (eigenvalues of \eqn{D - W}) and \code{car_eig}
#'   (eigenvalues of \eqn{D^{-1/2} W D^{-1/2}}).
#'
#' @keywords internal
prepare_sparse_car_W<-function(W){
  W<-as.matrix(W)
  # Check that W is a square symmetric adjacency matrix.
  if(nrow(W)!=ncol(W))stop("W must be a square matrix.")
  if(any(W!=t(W)))stop("W must be symmetric.")
  # The Stan CAR uses an unweighted graph with D_sparse = neighbour count.
  if(any(!W%in%c(0,1)))stop("W must be binary (0/1). Use W <- 1*(W>0).")
  # Remove any self-connections.
  diag(W)<-0
  # Number of pathway neighbors for each gene.
  degree<-rowSums(W)
  # Genes with at least one pathway neighbor are included in CAR.
  car_index<-which(degree>0)
  # Genes with no pathway neighbors receive an independent prior.
  iso_index<-which(degree==0)
  P_car<-length(car_index)
  P_iso<-length(iso_index)
  # Extract the connected portion of the pathway graph.
  W_car<-W[car_index,car_index,drop=FALSE]
  # Degree within the connected CAR graph.
  D_sparse<-rowSums(W_car)
  # Store each undirected edge only once.
  edge_index<-which(upper.tri(W_car)&W_car!=0,arr.ind=TRUE)
  W_n<-nrow(edge_index)
  W_sparse<-matrix(as.integer(edge_index),nrow=W_n,ncol=2)
  if(P_car>0){
    # Laplacian eigenvalues for the Leroux CAR log determinant:
    # log|(1-rho)I + rho(D - W)| = sum log(1 + rho(lap_eig - 1)).
    # Smallest is 0 per connected component (numerically ~ -1e-15): clamp at 0.
    lap_eig<-eigen(diag(D_sparse,nrow=P_car)-W_car,symmetric=TRUE,only.values=TRUE)$values
    lap_eig<-pmax(lap_eig,0)
    # Eigenvalues of D^{-1/2} W D^{-1/2} for the sparse proper CAR log determinant:
    # log|D - rho W| = sum log D + sum log(1 - rho car_eig).
    Dm12<-1/sqrt(D_sparse)
    car_eig<-eigen(W_car*outer(Dm12,Dm12),symmetric=TRUE,only.values=TRUE)$values
  }else{
    lap_eig<-numeric(0)
    car_eig<-numeric(0)
  }
  # Return everything required by Stan.
  list(
    P_car=P_car,
    car_index=as.integer(car_index),
    P_iso=P_iso,
    iso_index=as.integer(iso_index),
    W_n=W_n,
    W_sparse=W_sparse,
    D_sparse=as.numeric(D_sparse),
    lap_eig=as.numeric(lap_eig),
    car_eig=as.numeric(car_eig)
  )
}


## ------------------------------------------------------------
## Function to make knots, distances and RBF basis
## ------------------------------------------------------------

#' Construct Knots, Distances and RBF Basis
#'
#' Selects \code{r} knots by a seeded k-means and returns the spot-to-knot
#' distances and the RBF basis matrix built from the same knots, so \code{D}
#' and \code{Phi} always match.
#'
#' @param coords A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param r Integer. Number of knots; must be smaller than the number of
#'   spots (default: 30).
#' @param lengthscale Numeric or NULL. RBF lengthscale. If NULL, uses
#'   \code{ls_mult} times the median nearest-knot distance (default: NULL).
#' @param ls_mult Numeric. Multiplier for the default lengthscale
#'   (default: 1.5).
#' @param seed Integer or NULL. Seed for k-means (default: NULL).
#' @param nstart Integer. Number of k-means starts (default: 10).
#'
#' @return A named list with \code{knots}, \code{D} (\eqn{n \times r}),
#'   \code{Phi} (\eqn{n \times r}) and \code{lengthscale}.
#'
#' @importFrom stats kmeans dist median
#' @keywords internal
make_rbf_knots <- function(coords, r = 30, lengthscale = NULL, ls_mult = 1.5,
                           seed = NULL, nstart = 10) {
  coords <- as.matrix(coords)
  N <- nrow(coords)
  if (r >= N) stop("r must be smaller than the number of spots.")
  if (!is.null(seed)) set.seed(seed)
  km <- kmeans(coords, centers = r, nstart = nstart, iter.max = 100)
  knots <- km$centers
  D <- as.matrix(dist(rbind(coords, knots)))[1:N, (N + 1):(N + r), drop = FALSE]
  if (is.null(lengthscale)) {
    dk <- as.matrix(dist(knots))
    diag(dk) <- Inf
    lengthscale <- ls_mult * median(apply(dk, 1, min))
  }
  Phi <- exp(-(D^2) / (2 * lengthscale^2))
  list(knots = knots, D = D, Phi = Phi, lengthscale = lengthscale)
}


## ------------------------------------------------------------
## Function to check the pathway -> SE link
## ------------------------------------------------------------

#' Summarize the Pathway-SE Link in a Simulated Dataset
#'
#' Checks how strongly a simulated dataset or design ties SE status to the
#' pathway network.
#'
#' @param d A simulated dataset or design.
#' @param verbose Logical. If TRUE, prints the summary (default: TRUE).
#'
#' @return Invisibly, a named list of summaries including \code{n_SE},
#'   \code{SE_share}, neighbour SE shares, \code{nbr_corr},
#'   \code{cor_u_strength} and a per-pathway table.
#'
#' @importFrom stats cor
#' @keywords internal
summarize_pathway_se <- function(d, verbose = TRUE) {
  Z <- d$Z
  W <- d$W
  deg <- rowSums(W)
  M <- d$pathway_membership
  has_nb <- deg > 0
  nb_se <- rep(NA_real_, length(Z))
  nb_se[has_nb] <- as.vector(W %*% Z)[has_nb] / deg[has_nb]
  pattern <- if (!is.null(d$eta)) apply(d$eta, 2, sd) > 0 else d$sigma_eta > 0
  se <- Z == 1
  per_pathway <- data.frame(
    pathway = colnames(M), n_genes = colSums(M), n_SE = colSums(M * Z),
    SE_share = colSums(M * Z) / colSums(M), row.names = NULL
  )
  u_ok <- !is.null(d$u_centered) && all(is.finite(d$u_centered))
  out <- list(
    n_SE = sum(Z), SE_share = mean(Z),
    SE_with_pattern = sum(pattern[se]), nonSE_with_pattern = sum(pattern[!se]),
    nb_SE_share_SEgenes = mean(nb_se[se], na.rm = TRUE),
    nb_SE_share_nonSEgenes = mean(nb_se[!se], na.rm = TRUE),
    pathways_with_SE = sum(per_pathway$n_SE > 0),
    nbr_corr = if (is.null(d$nbr_corr)) NA_real_ else d$nbr_corr,
    cor_u_strength = if (u_ok && sum(se) > 2) cor(d$u_centered[se], d$sigma_eta[se]) else NA_real_,
    per_pathway = per_pathway
  )
  if (verbose) {
    cat("SE genes:", out$n_SE, sprintf("(%.0f%%)", 100 * out$SE_share), "\n")
    cat(
      "Genes with a spatial pattern: SE", out$SE_with_pattern, "of", out$n_SE,
      "| non-SE", out$nonSE_with_pattern, "\n"
    )
    cat(sprintf(
      "Share of neighbours that are SE: SE genes %.2f vs non-SE genes %.2f\n",
      out$nb_SE_share_SEgenes, out$nb_SE_share_nonSEgenes
    ))
    cat("Pathways containing SE genes:", out$pathways_with_SE, "of", nrow(per_pathway), "\n")
    cat("SE share per pathway:", paste(sprintf("%.0f%%", 100 * sort(per_pathway$SE_share)), collapse = " "), "\n")
    if (u_ok) cat(sprintf("corr(pathway effect u, pattern strength) among SE genes: %.2f\n", out$cor_u_strength))
  }
  invisible(out)
}


## ------------------------------------------------------------
## Function to validate exclusive pathways
## ------------------------------------------------------------

#' Validate an Exclusive Pathway Structure
#'
#' Verifies that SE and non-SE genes belong exclusively to their designated
#' pathways.
#'
#' @param M Binary gene-by-pathway membership matrix.
#' @param Z Integer vector of SE indicators (1 = SE).
#' @param se_pathways Integer or character vector. SE pathways.
#'
#' @return A named list with \code{valid}, \code{checks} and a per-pathway
#'   \code{summary}.
#'
#' @keywords internal
validate_exclusive_pathways <- function(M, Z, se_pathways) {
  M <- as.matrix(M)
  P <- nrow(M)
  K <- ncol(M)
  if (length(Z) != P) stop("Z must have one value per gene.")
  if (is.character(se_pathways)) se_pathways <- match(se_pathways, colnames(M))
  se_pathways <- sort(unique(as.integer(se_pathways)))
  if (length(se_pathways) == 0 || anyNA(se_pathways) || any(se_pathways < 1 | se_pathways > K)) stop("Invalid se_pathways.")
  nonse_pathways <- setdiff(seq_len(K), se_pathways)
  se_idx <- which(Z == 1)
  nonse_idx <- which(Z == 0)
  checks <- c(
    every_gene_has_pathway = all(rowSums(M) > 0),
    every_pathway_has_gene = all(colSums(M) > 0),
    se_only_in_se_pathways = all(M[se_idx, nonse_pathways, drop = FALSE] == 0),
    nonse_only_in_nonse_pathways = all(M[nonse_idx, se_pathways, drop = FALSE] == 0),
    all_se_pathways_have_se = all(colSums(M[se_idx, se_pathways, drop = FALSE]) > 0),
    all_nonse_pathways_have_nonse = all(colSums(M[nonse_idx, nonse_pathways, drop = FALSE]) > 0)
  )
  summary <- data.frame(
    pathway = colnames(M),
    type = ifelse(seq_len(K) %in% se_pathways, "SE", "NonSE"),
    total_genes = colSums(M),
    se_genes = colSums(M[se_idx, , drop = FALSE]),
    nonse_genes = colSums(M[nonse_idx, , drop = FALSE])
  )
  list(valid = all(checks), checks = checks, summary = summary)
}


## ------------------------------------------------------------
## Function to validate exclupath pathways
## ------------------------------------------------------------

#' Validate an Exclupath Pathway Structure
#'
#' Checks that SE genes appear only in SE pathways, every gene and pathway is
#' used, and every SE pathway contains at least one non-SE gene.
#'
#' @param M Binary gene-by-pathway membership matrix.
#' @param Z Integer vector of SE indicators (1 = SE).
#' @param se_pathways Integer vector. SE pathways.
#'
#' @return A named list with \code{valid} and \code{checks}.
#'
#' @keywords internal
validate_exclupath <- function(M, Z, se_pathways) {
  M <- as.matrix(M); se_pathways <- as.integer(se_pathways); non_pw <- setdiff(seq_len(ncol(M)), se_pathways)
  checks <- c(every_gene_has_pathway = all(rowSums(M) > 0),
              every_pathway_has_gene = all(colSums(M) > 0),
              se_only_in_se_pathways = all(M[Z == 1, non_pw, drop = FALSE] == 0),
              nonse_pathways_have_no_se = all(colSums(M[Z == 1, non_pw, drop = FALSE]) == 0),
              # [CHANGED 5] was "SE pathways mostly SE"; now "at least one non-SE gene in every SE pathway"
              se_pathways_have_nonse = all(colSums(M[Z == 0, se_pathways, drop = FALSE]) >= 1))
  list(valid = all(checks), checks = checks)
}


## ------------------------------------------------------------
## Function to summarize pathway composition
## ------------------------------------------------------------

#' Pathway Composition
#'
#' One row per pathway: number of genes, SE genes, non-SE genes, percent SE
#' and its type. Works for every design and on saved replicates.
#'
#' @param M Binary gene-by-pathway membership matrix.
#' @param Z Integer vector of SE indicators (1 = SE).
#' @param se_pathways Integer or character vector, or NULL. Designated SE
#'   pathways, if known (default: NULL).
#'
#' @return A data frame with one row per pathway. \code{attr(., "genes")}
#'   counts genes by number of pathways.
#'
#' @keywords internal
pathway_composition <- function(M, Z, se_pathways = NULL) {
  M <- as.matrix(M); Z <- as.integer(Z)
  if (is.null(colnames(M))) colnames(M) <- paste0("Pathway", seq_len(ncol(M)))
  out <- data.frame(pathway = colnames(M), n_genes = colSums(M), n_SE = colSums(M[Z == 1, , drop = FALSE]),
                    n_nonSE = colSums(M[Z == 0, , drop = FALSE]), row.names = NULL)
  out$pct_SE <- round(100 * out$n_SE / pmax(out$n_genes, 1), 1)
  out$type <- if (!is.null(se_pathways) && length(se_pathways)) {
    sp <- if (is.character(se_pathways)) match(se_pathways, colnames(M)) else as.integer(se_pathways)
    ifelse(seq_len(ncol(M)) %in% sp, "SE pathway", "non-SE pathway")
  } else ifelse(out$n_SE == out$n_genes, "all SE", ifelse(out$n_SE == 0, "no SE", "mixed"))
  out <- out[order(-out$pct_SE, out$pathway), ]
  k <- rowSums(M)
  attr(out, "genes") <- as.data.frame.matrix(table(status = ifelse(Z == 1, "SE", "non-SE"),
                                                   pathways_per_gene = factor(k, levels = sort(unique(k)))))
  out
}


## ------------------------------------------------------------
## Function to print pathway composition
## ------------------------------------------------------------

#' Print Pathway Composition
#'
#' Prints the output of \code{pathway_composition}.
#'
#' @param pc Output of \code{pathway_composition}.
#'
#' @return Invisibly, \code{pc}.
#'
#' @keywords internal
print_composition <- function(pc) {
  cat("\n==== PATHWAY COMPOSITION ====\n"); print(pc, row.names = FALSE)
  cat("\nGenes by number of pathways they belong to:\n"); print(attr(pc, "genes"))
  cat(sprintf("\nTotal: %d SE and %d non-SE genes\n", sum(attr(pc, "genes")["SE", ]), sum(attr(pc, "genes")["non-SE", ])))
  invisible(pc)
}


## ------------------------------------------------------------
## Function to label a scenario for file names
## ------------------------------------------------------------

#' Scenario Label for File Names
#'
#' Short label of a simulation scenario, e.g. \code{"car-rho0.90-su1.5"},
#' \code{"hot2-nse20"} or \code{"random-nse20"}. Pattern-based designs get
#' \code{"-pattern"} plus any non-default pattern setting.
#'
#' @param param The \code{param} list of a design or dataset.
#'
#' @return A character string.
#'
#' @keywords internal
design_label <- function(param) {
  se_design <- if (is.null(param$se_design)) "car" else param$se_design
  base <- switch(se_design,
                 "car" = sprintf("car-rho%s-su%s", formatC(param$rho_u, format = "f", digits = 2), format(param$sigma_u)),
                 "hot_pathways" = paste0(
                   sprintf("hot%d-nse%d", param$n_hot, param$n_se),
                   if (isTRUE(param$prop_outside > 0)) sprintf("-out%02d", round(100 * param$prop_outside)) else ""
                 ),
                 "random" = sprintf("random-nse%d", param$n_se),
                 "exclusive_pathways" = sprintf("exclusive%d-nse%d", param$n_se_pathways, param$n_se),
                 "exclupath" = sprintf("exclupath%d-nse%d-se%g-ns%g", param$n_se_pathways, param$n_se, param$se_per_pathway, param$nonse_per_se_pathway)
  )
  if (identical(param$signal_model, "pattern")) {
    pt <- param$pattern_types
    base <- paste0(base, "-pattern",
                   if (!is.null(pt) && !identical(pt, "mixed")) paste0("-", paste(pt, collapse = "+")) else "",
                   if (identical(param$pathway_patterns, "shared")) "-shared" else "",
                   if (!is.null(param$signal_scenario) && param$signal_scenario != "mixed") paste0("-", param$signal_scenario) else "")
  }
  base
}


## ------------------------------------------------------------
## Function to compute rank-based AUC
## ------------------------------------------------------------

#' Rank-Based AUC
#'
#' @param score Numeric vector. Detection score (higher = more likely SE).
#' @param truth Integer vector. True SE indicators (1 = SE).
#'
#' @return Numeric AUC, or NA if only one class is present.
#'
#' @keywords internal
auc_rank <- function(score, truth) {
  truth <- as.integer(truth)
  n1 <- sum(truth)
  n0 <- sum(1 - truth)
  if (n1 == 0 || n0 == 0) {return(NA_real_)}
  r <- rank(score)
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)}


## ------------------------------------------------------------
## Function to compute average precision
## ------------------------------------------------------------

#' Average Precision
#'
#' Area under the precision-recall curve.
#'
#' @param score Numeric vector. Detection score (higher = more likely SE).
#' @param truth Integer vector. True SE indicators (1 = SE).
#'
#' @return Numeric average precision, or NA if there are no SE genes.
#'
#' @keywords internal
avg_precision <- function(score, truth) {
  truth <- as.integer(truth)
  if (sum(truth) == 0) {return(NA_real_)}
  o <- order(score, decreasing = TRUE)
  t <- truth[o]
  sum((cumsum(t) / seq_along(t))[t == 1]) / sum(t)
}


## ------------------------------------------------------------
## Function to select genes by Bayesian FDR
## ------------------------------------------------------------

#' Bayesian FDR Selection
#'
#' Selects the largest set of genes whose expected false discovery rate,
#' based on posterior SE probabilities, is at most \code{alpha}.
#'
#' @param prob Numeric vector. Posterior SE probabilities.
#' @param alpha Numeric. Target Bayesian FDR (default: 0.05).
#'
#' @return Integer vector of calls (1 = SE).
#'
#' @keywords internal
bayes_fdr_select <- function(prob, alpha = 0.05) {
  o <- order(prob, decreasing = TRUE)
  efdr <- cumsum(1 - prob[o]) / seq_along(o)
  k <- sum(efdr <= alpha)
  sel <- integer(length(prob))
  if (k > 0) sel[o[seq_len(k)]] <- 1L
  sel}


## ------------------------------------------------------------
## Function to compute Matthews correlation coefficient
## ------------------------------------------------------------

#' Matthews Correlation Coefficient
#'
#' @param TP,FP,TN,FN Integer. Confusion-matrix counts.
#'
#' @return Numeric MCC (0 if undefined).
#'
#' @keywords internal
mcc_score <- function(TP, FP, TN, FN) {
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (den == 0) return(0)
  ((TP * TN) - (FP * FN)) / den
}


## ------------------------------------------------------------
## Function to compute SPHERE detection metrics
## ------------------------------------------------------------

#' SPHERE Detection Metrics
#'
#' Computes AUC, average precision and, at both a probability threshold and a
#' Bayesian FDR cut-off, TPR, FPR, precision, FDR, F1 and MCC. If
#' \code{signal_class} is given, also reports TPR by weak / moderate / strong
#' signal.
#'
#' @param prob Numeric vector. Posterior SE probabilities.
#' @param truth Integer vector. True SE indicators (1 = SE).
#' @param signal_class Character vector or NULL. Signal class per gene
#'   (default: NULL).
#' @param threshold Numeric. Probability threshold for calls (default: 0.5).
#' @param fdr_alpha Numeric. Target Bayesian FDR (default: 0.05).
#'
#' @return A one-row data frame of metrics.
#'
#' @importFrom stats setNames
#' @keywords internal
detection_metrics <- function(prob, truth, signal_class = NULL, threshold = 0.5, fdr_alpha = 0.05) {
  truth <- as.integer(truth)
  rates <- function(pred) {
    TP <- sum(pred == 1 & truth == 1)
    FP <- sum(pred == 1 & truth == 0)
    TN <- sum(pred == 0 & truth == 0)
    FN <- sum(pred == 0 & truth == 1)
    TPR <- TP / (TP + FN)
    FPR <- FP / (FP + TN)
    prec <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
    # Same conventions as evaluate_method(): nothing called -> FDR 0, F1 0, Precision NA.
    F1 <- 2 * TP / (2 * TP + FP + FN)
    c(TPR = TPR, FPR = FPR, Precision = prec,
      FDR = if (TP + FP > 0) FP / (TP + FP) else 0,
      F1 = F1, MCC = mcc_score(TP, FP, TN, FN))
  }
  thr <- rates(as.integer(prob >= threshold))
  sel <- bayes_fdr_select(prob, fdr_alpha)
  bf <- rates(sel)
  out <- data.frame(
    AUC = auc_rank(prob, truth), AP = avg_precision(prob, truth),
    t(setNames(thr, paste0(names(thr), "_thr"))),
    t(setNames(bf, paste0(names(bf), "_bfdr"))), check.names = FALSE
  )
  if (!is.null(signal_class)) {
    for (cl in c("weak", "moderate", "strong")) {
      idx <- which(signal_class == cl)
      out[[paste0("TPR_bfdr_", cl)]] <- if (length(idx) > 0) mean(sel[idx] == 1) else NA_real_
    }
  }
  out
}
