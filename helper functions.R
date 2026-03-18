##----------------------------
## Plot all Figures
##----------------------------

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










