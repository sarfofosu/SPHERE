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
## Function to filter genes with less total counts
## ------------------------------------------------------------

#' Filter Genes with Low Total Counts
#'
#' Removes genes (columns) whose total count across all spots falls
#' below a specified threshold.
#'
#' @param df A matrix or data frame of gene expression counts.
#' @param threshold Numeric. Minimum total count required to keep a
#'   gene (default: 10).
#' @param keep_first Logical. If TRUE, always keeps the first column
#'   regardless of its sum (default: TRUE).
#'
#' @return A filtered \code{data.table} with low-count genes removed.
#'
#' @importFrom data.table as.data.table
#' @export
filter_columns_by_sum <- function(df, threshold = 10, keep_first = TRUE) {
  ..cols_keep <- NULL  # fix R CMD CHECK note
  dt <- as.data.table(df)
  if (keep_first) {
    numeric_part <- dt[, -1, with = FALSE]
    sums         <- colSums(numeric_part, na.rm = TRUE)
    cols_keep    <- c(names(dt)[1], names(numeric_part)[sums > threshold])
  } else {
    sums      <- colSums(dt, na.rm = TRUE)
    cols_keep <- names(dt)[sums > threshold]
  }
  return(dt[, ..cols_keep])
}


## ------------------------------------------------------------
## Function to make RBF basis matrix
## ------------------------------------------------------------

#' Construct Low-Rank RBF Gaussian Process Basis Matrix
#'
#' Computes a radial basis function (RBF) basis matrix using \code{r}
#' knot centers selected by k-means clustering. Used to construct the
#' low-rank GP approximation in SPHERE.
#'
#' @param coords A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param r Integer. Number of knots (inducing points) (default: 30).
#' @param lengthscale Numeric or NULL. RBF lengthscale. If NULL, uses
#'   the median pairwise distance between knots (default: 1.5).
#'
#' @return A numeric matrix of dimension \eqn{n \times r} containing
#'   the RBF basis function values.
#'
#' @importFrom stats kmeans dist median
#' @export
make_rbf_basis <- function(coords, r = 30, lengthscale = 1.5) {
  N  <- nrow(coords)
  km <- kmeans(coords, centers = r)
  knots <- km$centers
  d  <- as.matrix(dist(rbind(coords, knots)))
  D  <- d[1:N, (N + 1):(N + r)]
  if (is.null(lengthscale)) {
    lengthscale <- median(as.matrix(dist(knots)))
  }
  Phi <- exp(-(D^2) / (2 * lengthscale^2))
  return(Phi)
}


## ------------------------------------------------------------
## Function to make RBF distance matrix
## ------------------------------------------------------------

#' Compute Spot-to-Knot Distance Matrix for GP Approximation
#'
#' Computes the matrix of Euclidean distances from each spatial spot
#' to each of \code{r} knot centers selected by k-means clustering.
#'
#' @param coords A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param r Integer. Number of knots (inducing points) (default: 30).
#'
#' @return A numeric matrix of dimension \eqn{n \times r} containing
#'   Euclidean distances from spots to knots.
#'
#' @importFrom stats kmeans dist
#' @export
make_rbf_dist <- function(coords, r = 30) {
  coords <- as.matrix(coords)
  N      <- nrow(coords)
  km     <- kmeans(coords, centers = r)
  knots  <- km$centers
  all_d  <- as.matrix(dist(rbind(coords, knots)))
  D      <- all_d[1:N, (N + 1):(N + r)]
  D
}


## ------------------------------------------------------------
## Function to generate gene pathway group
## ------------------------------------------------------------

#' Generate Gene Pathway Group Membership
#'
#' Assigns genes to \code{G} pathway groups in a round-robin fashion.
#' Genes are assigned sequentially: gene 1 to group 1, gene 2 to group
#' 2, ..., gene G to group G, gene G+1 to group 1, and so on.
#'
#' @param p Integer. Total number of genes.
#' @param G Integer. Number of pathway groups (default: 3).
#'
#' @return An integer vector of length \code{p} giving pathway group
#'   membership for each gene.
#'
#' @export
generate_gene_grp <- function(p, G = 3) {
  reps <- ceiling(p / G)
  rep(1:G, each = reps)[1:p]
}


## ------------------------------------------------------------
## Function to get the tau and ell for SE genes
## ------------------------------------------------------------

#' Get Tau and Ell Parameters for Spatially Expressed Genes
#'
#' Extracts the GP amplitude (tau) and lengthscale (ell) parameters
#' for spatially expressed (SE) genes based on a given proportion.
#'
#' @param p Integer. Total number of genes.
#' @param prop Numeric vector of length 2 giving the proportion of
#'   non-SE and SE genes respectively (e.g. \code{c(0.7, 0.3)}).
#' @param t_gs Numeric vector of GP amplitude values for SE genes.
#' @param l_gs Numeric vector of GP lengthscale values for SE genes.
#'
#' @return A named list with elements \code{tau_gs} and \code{ell_gs},
#'   each a numeric vector of length \code{round(prop[2] * p)}.
#'
#' @export
get_tau_ell <- function(p, prop, t_gs, l_gs) {
  se_count <- round(prop[2] * p)
  list(
    tau_gs = t_gs[1:se_count],
    ell_gs = l_gs[1:se_count]
  )
}


## ------------------------------------------------------------
## Relative expression
## ------------------------------------------------------------

#' Compute Relative Gene Expression
#'
#' Applies a log1p transformation followed by min-max normalization
#' to convert raw counts to relative expression values in \eqn{[0, 1]}.
#'
#' @param raw_exp A numeric vector or matrix of raw expression counts.
#'
#' @return A numeric vector or matrix of relative expression values
#'   scaled to \eqn{[0, 1]}.
#'
#' @export
relative_expr <- function(raw_exp) {
  raw_exp_log <- log1p(raw_exp)
  (raw_exp_log - min(raw_exp_log)) / (max(raw_exp_log) - min(raw_exp_log))
}


## ------------------------------------------------------------
## Create pathway effects (internal)
## ------------------------------------------------------------

#' Simulate CAR Pathway Effects for Genes
#'
#' Simulates gene-level pathway effects \eqn{\beta_j} from a
#' Conditional Autoregressive (CAR) prior. Genes in the same pathway
#' group are treated as spatial neighbors and share information.
#'
#' @param num_genes Integer. Total number of genes.
#' @param gene_grp Integer vector of length \code{num_genes} giving
#'   pathway group membership for each gene.
#' @param rho Numeric. CAR spatial correlation parameter (default: 0.9).
#' @param tau_beta Numeric. CAR precision parameter (default: 1).
#'
#' @return A numeric vector of length \code{num_genes} containing
#'   simulated gene-level pathway effects.
#'
#' @importFrom stats rnorm
#' @keywords internal
create_car_beta2 <- function(num_genes, gene_grp, rho = 0.9, tau_beta = 1) {
  stopifnot(length(gene_grp) == num_genes)

  W <- matrix(0, num_genes, num_genes)
  for (i in 1:(num_genes - 1)) {
    for (j in (i + 1):num_genes) {
      if (gene_grp[i] == gene_grp[j]) {
        W[i, j] <- 1
        W[j, i] <- 1
      }
    }
  }

  deg  <- rowSums(W)
  D_mat <- diag(deg)
  Q0   <- D_mat - rho * W

  iso <- which(deg == 0)
  if (length(iso) > 0) Q0[iso, iso] <- 1

  Q  <- tau_beta * (Q0 + t(Q0)) / 2
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0)
    stop(sprintf("Q is not positive definite; min eigenvalue = %.3e. Reduce rho.",
                 min(ev)))

  R    <- chol(Q)
  z    <- rnorm(num_genes)
  Beta <- backsolve(R, z, transpose = TRUE)
  Beta <- backsolve(R, Beta)
  Beta
}


## ------------------------------------------------------------
## Create gene-specific spatial variability (internal)
## ------------------------------------------------------------

#' Simulate Gene-Specific Gaussian Process Spatial Effects.
#'
#' Simulates spatial expression effects for SE genes using
#' gene-specific Gaussian Process covariance functions with an
#' RBF (squared exponential) kernel.
#'
#' @param spots A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param num_genes Integer. Number of SE genes to simulate.
#' @param tau_gs Numeric vector of length \code{num_genes}. GP amplitude
#'   (signal variance) for each gene.
#' @param ell_gs Numeric vector of length \code{num_genes}. GP lengthscale
#'   for each gene.
#' @param nugget Numeric. Nugget added to diagonal for numerical stability
#'   (default: 1e-6).
#'
#' @return A numeric matrix of dimension \eqn{n \times} \code{num_genes}
#'   containing simulated GP spatial effects.
#'
#' @importFrom stats dist rnorm
#' @keywords internal
sim_gs_eta2 <- function(spots, num_genes, tau_gs, ell_gs, nugget = 1e-6) {
  n <- nrow(spots)
  if (length(tau_gs) != num_genes || length(ell_gs) != num_genes)
    stop("tau_gs and ell_gs must have length num_genes")

  c_mat   <- matrix(0, n, num_genes)
  dist_sq <- as.matrix(dist(spots))^2

  for (j in 1:num_genes) {
    K_j        <- tau_gs[j] * exp(-dist_sq / (2 * ell_gs[j]^2))
    diag(K_j)  <- diag(K_j) + nugget
    R          <- chol(K_j)
    z          <- rnorm(n)
    c_mat[, j] <- t(R) %*% z
  }
  return(c_mat)
}
