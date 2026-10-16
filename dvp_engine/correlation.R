# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - correlation
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/correlations/helpers/corr_helpers.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

compute_correlation <- function(data, 
                                 method = "pearson", 
                                 use = "pairwise.complete.obs",
                                 p_adjust = NULL) {
  
  # Convert to matrix if needed
  if (is.data.frame(data)) {
    data <- as.matrix(data[, sapply(data, is.numeric)])
  }
  
  n_vars <- ncol(data)
  var_names <- colnames(data)
  
  # Compute correlation matrix
  corr_mat <- cor(data, method = method, use = use)
  
  # Initialize p-value and sample size matrices
  p_mat <- matrix(NA, nrow = n_vars, ncol = n_vars,
                  dimnames = list(var_names, var_names))
  n_mat <- matrix(NA, nrow = n_vars, ncol = n_vars,
                  dimnames = list(var_names, var_names))
  
  # Compute p-values and sample sizes
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      # Get complete cases for this pair
      complete_idx <- complete.cases(data[, c(i, j)])
      n_complete <- sum(complete_idx)
      n_mat[i, j] <- n_complete
      n_mat[j, i] <- n_complete
      
      if (n_complete >= 3) {
        test <- cor.test(data[complete_idx, i], data[complete_idx, j], method = method)
        p_mat[i, j] <- test$p.value
        p_mat[j, i] <- test$p.value
      }
    }
  }
  
  # Set diagonal values
  diag(p_mat) <- 0
  diag(n_mat) <- apply(data, 2, function(x) sum(!is.na(x)))
  
  # Adjust p-values if requested
  if (!is.null(p_adjust)) {
    p_upper <- p_mat[upper.tri(p_mat)]
    p_adjusted <- p.adjust(p_upper, method = p_adjust)
    p_mat_adjusted <- p_mat
    p_mat_adjusted[upper.tri(p_mat_adjusted)] <- p_adjusted
    p_mat_adjusted[lower.tri(p_mat_adjusted)] <- t(p_mat_adjusted)[lower.tri(p_mat_adjusted)]
    p_mat <- p_mat_adjusted
  }
  
  list(
    correlation = corr_mat,
    p_value = p_mat,
    n = n_mat,
    method = method,
    n_vars = n_vars,
    var_names = var_names
  )
}

