# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - pca_helpers
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/pcabiplot/helpers/pca_helpers.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

standardize_prcomp <- function(pca) {
  
  eigenvalues <- pca$sdev^2
  total_var <- sum(eigenvalues)
  variance_percent <- (eigenvalues / total_var) * 100
  cumulative_variance <- cumsum(variance_percent)
  
  # Scores (individuals)
  scores <- pca$x
  rownames(scores) <- rownames(scores) %||% paste0("Ind", 1:nrow(scores))
  
  # Loadings (variables)
  loadings <- pca$rotation
  rownames(loadings) <- rownames(loadings) %||% paste0("Var", 1:nrow(loadings))
  
  # Calculate contributions (approximation)
  contrib_var <- sweep(loadings^2, 2, colSums(loadings^2), "/") * 100
  contrib_ind <- sweep(scores^2, 2, colSums(scores^2), "/") * 100
  
  # Calculate cos2 (quality of representation)
  row_norms_var <- rowSums(loadings^2)
  cos2_var <- sweep(loadings^2, 1, row_norms_var, "/")
  
  row_norms_ind <- rowSums(scores^2)
  cos2_ind <- sweep(scores^2, 1, row_norms_ind, "/")
  
  list(
    type = "prcomp",
    eigenvalues = eigenvalues,
    variance_percent = variance_percent,
    cumulative_variance = cumulative_variance,
    scores = scores,
    loadings = loadings,
    contributions_var = contrib_var,
    contributions_ind = contrib_ind,
    cos2_var = cos2_var,
    cos2_ind = cos2_ind,
    center = pca$center,
    scale = pca$scale,
    n_components = ncol(scores),
    n_individuals = nrow(scores),
    n_variables = nrow(loadings)
  )
}

standardize_factominer <- function(pca) {
  
  eigenvalues <- pca$eig[, 1]
  variance_percent <- pca$eig[, 2]
  cumulative_variance <- pca$eig[, 3]
  
  # Scores (individuals)
  scores <- pca$ind$coord
  
  # Loadings (variables)
  loadings <- pca$var$coord
  
  # Contributions
  contrib_var <- pca$var$contrib
  contrib_ind <- pca$ind$contrib
  
  # Cos2
  cos2_var <- pca$var$cos2
  cos2_ind <- pca$ind$cos2
  
  list(
    type = "factominer",
    eigenvalues = eigenvalues,
    variance_percent = variance_percent,
    cumulative_variance = cumulative_variance,
    scores = scores,
    loadings = loadings,
    correlations = pca$var$cor %||% loadings,   # PCA-002: true variable-PC correlations
    contributions_var = contrib_var,
    contributions_ind = contrib_ind,
    cos2_var = cos2_var,
    cos2_ind = cos2_ind,
    quali_sup = pca$quali.sup,
    quanti_sup = pca$quanti.sup,
    n_components = ncol(scores),
    n_individuals = nrow(scores),
    n_variables = nrow(loadings)
  )
}

