# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - pca
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/pcabiplot/server/pca_computation.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

compute_pca <- function(data,
                         n_components = 5,
                         scale = TRUE,
                         center = TRUE,
                         quali_sup = NULL,
                         quanti_sup = NULL,
                         ind_sup = NULL,
                         method = "factominer") {
  
  # Ensure n_components doesn't exceed dimensions
  max_components <- min(nrow(data) - 1, ncol(data))
  n_components <- min(n_components, max_components)
  
  if (method == "factominer" && requireNamespace("FactoMineR", quietly = TRUE)) {
    
    # Use FactoMineR
    pca_args <- list(
      X = data,
      ncp = n_components,
      scale.unit = scale,
      graph = FALSE
    )
    
    # Add supplementary variables if specified
    if (!is.null(quali_sup) && length(quali_sup) > 0) {
      pca_args$quali.sup <- quali_sup
    }
    if (!is.null(quanti_sup) && length(quanti_sup) > 0) {
      pca_args$quanti.sup <- quanti_sup
    }
    if (!is.null(ind_sup) && length(ind_sup) > 0) {
      pca_args$ind.sup <- ind_sup
    }
    
    pca_result <- do.call(FactoMineR::PCA, pca_args)
    
    # Standardize result
    pca_standard <- standardize_factominer(pca_result)
    
  } else {
    
    # Use base R prcomp
    pca_result <- prcomp(
      data,
      scale. = scale,
      center = center,
      rank. = n_components
    )
    
    # Standardize result
    pca_standard <- standardize_prcomp(pca_result)
    
  }
  
  # Add metadata
  pca_standard$n_components <- n_components
  pca_standard$method <- method
  pca_standard$scale <- scale
  pca_standard$center <- center
  pca_standard$n_obs <- nrow(data)
  pca_standard$n_vars <- ncol(data)
  pca_standard$var_names <- colnames(data)
  pca_standard$obs_names <- rownames(data) %||% paste0("Obs", 1:nrow(data))
  
  return(pca_standard)
}

