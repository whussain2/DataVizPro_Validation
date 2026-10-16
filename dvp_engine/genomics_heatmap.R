# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - genomics_heatmap
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/genomicvisuals/genomics_heatmaps/helpers/gh_helpers.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

GH_CHUNK_COLS <- 4000L

.gh_col_mode <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  tb <- table(v)
  as.numeric(names(tb)[which.max(tb)])
}

gh_impute <- function(M, method = "mean") {
  M <- as.matrix(M)
  if (!anyNA(M) || identical(method, "keep")) return(M)
  if (identical(method, "complete")) {
    return(M[stats::complete.cases(M), , drop = FALSE])
  }
  na <- is.na(M)
  fill <- switch(method,
    mode   = apply(M, 2, .gh_col_mode),
    median = apply(M, 2, stats::median, na.rm = TRUE),
    zero   = rep(0, ncol(M)),
    colMeans(M, na.rm = TRUE))                      # mean (default)
  fill[!is.finite(fill)] <- 0
  # Column-wise fill through the NA index (no per-column loop on big matrices)
  idx <- which(na, arr.ind = TRUE)
  M[idx] <- fill[idx[, 2]]
  M
}

.gh_chunk_accumulate <- function(M, fn, chunk = GH_CHUNK_COLS, progress = NULL) {
  n <- nrow(M); m <- ncol(M)
  acc <- matrix(0, n, n)
  starts <- seq(1L, m, by = chunk)
  for (s in starts) {
    e <- min(s + chunk - 1L, m)
    block <- fn(M[, s:e, drop = FALSE], s:e)
    acc[] <- acc + block
    if (is.function(progress)) progress(e / m)
  }
  acc
}

.gh_row_crossprod <- function(M, center = FALSE, progress = NULL) {
  if (!isTRUE(center))
    return(.gh_chunk_accumulate(M, function(X, cols) tcrossprod(X), progress = progress))

  n <- nrow(M); m <- ncol(M)
  origins <- M[, 1L]
  residual_sums <- numeric(n)
  starts <- seq(1L, m, by = GH_CHUNK_COLS)
  for (s in starts) {
    e <- min(s + GH_CHUNK_COLS - 1L, m)
    residual_sums <- residual_sums + rowSums(sweep(M[, s:e, drop = FALSE], 1, origins, "-"))
    if (is.function(progress)) progress(0.2 * e / m)
  }
  residual_means <- residual_sums / m
  .gh_chunk_accumulate(M, function(X, cols) {
    X <- sweep(X, 1, origins, "-")
    X <- sweep(X, 1, residual_means, "-")
    tcrossprod(X)
  }, progress = if (is.function(progress)) function(p) progress(0.2 + 0.8 * p) else NULL)
}

gh_kinship <- function(M, method = "vanraden", progress = NULL) {
  M <- as.matrix(M)
  if (!is.numeric(M)) storage.mode(M) <- "double"
  if (anyNA(M)) M <- gh_impute(M, "mean")
  n <- nrow(M)
  ids <- rownames(M) %||% paste0("A", seq_len(n))

  # Correlation / cosine / covariance work on ANY assay and are routed here so
  # the server has a single similarity entry point.
  if (method %in% c("pearson_sim", "spearman_sim", "cosine", "covariance", "centred_sim")) {
    S <- gh_similarity_matrix(M, method, progress = progress)
    dimnames(S) <- list(ids, ids)
    return(S)
  }

  if (identical(method, "ibs")) {
    K <- gh_ibs_similarity(M, progress = progress)
    dimnames(K) <- list(ids, ids)
    return(K)
  }

  p <- colMeans(M) / 2
  keep <- is.finite(p) & p > 0 & p < 1          # drop monomorphic markers
  if (!any(keep)) stop("All markers are monomorphic; kinship is undefined.")
  M <- M[, keep, drop = FALSE]; p <- p[keep]
  m <- ncol(M)

  if (identical(method, "vanraden")) {
    denom <- 2 * sum(p * (1 - p))
    K <- .gh_chunk_accumulate(M, function(X, cols) {
      Z <- sweep(X, 2, 2 * p[cols], "-")
      tcrossprod(Z)
    }, progress = progress) / denom
  } else {
    # yang / astle : standardise each marker by sqrt(2p(1-p))
    K <- .gh_chunk_accumulate(M, function(X, cols) {
      pc <- p[cols]
      W  <- sweep(sweep(X, 2, 2 * pc, "-"), 2, sqrt(2 * pc * (1 - pc)), "/")
      tcrossprod(W)
    }, progress = progress) / m
    if (identical(method, "astle")) K <- K / 2
  }
  dimnames(K) <- list(ids, ids)
  K
}

gh_similarity_matrix <- function(M, method = "pearson_sim", progress = NULL) {
  M <- as.matrix(M)
  if (anyNA(M)) M <- gh_impute(M, "mean")
  S <- switch(method,
    pearson_sim  = {
      cross <- .gh_row_crossprod(M, center = TRUE, progress = progress)
      scale <- sqrt(pmax(diag(cross), 0))
      cross / outer(scale, scale)
    },
    spearman_sim = stats::cor(t(M), method = "spearman"),
    centred_sim  = {
      cross <- .gh_row_crossprod(M, center = TRUE, progress = progress)
      scale <- sqrt(pmax(diag(cross), 0))
      cross / outer(scale, scale)
    },
    covariance   = .gh_row_crossprod(M, center = TRUE, progress = progress) /
                   max(1, ncol(M) - 1),
    cosine       = {
      nrm <- sqrt(rowSums(M^2)); nrm[!is.finite(nrm) | nrm == 0] <- 1
      out <- tcrossprod(M / nrm)
      if (is.function(progress)) progress(1)
      out
    },
    stats::cor(t(M)))
  S[!is.finite(S)] <- 0
  dimnames(S) <- list(rownames(M), rownames(M))
  S
}

gh_ibs_similarity <- function(M, progress = NULL) {
  M <- as.matrix(M)
  if (anyNA(M)) M <- gh_impute(M, "mean")
  M <- round(M); M[M < 0] <- 0; M[M > 2] <- 2
  m <- ncol(M)
  # sum|a-b| = sum (a-b)^2 - 2 * #pairs{|a-b| = 2}
  #          = (sq_i + sq_j - 2 M M') - 2 (A0 A2' + A2 A0')
  sq <- rowSums(M^2)
  S  <- .gh_chunk_accumulate(M, function(X, cols) {
    A0 <- (X == 0) * 1; A2 <- (X == 2) * 1
    C02 <- tcrossprod(A0, A2)
    -2 * tcrossprod(X) - 2 * (C02 + t(C02))
  }, progress = progress)
  S <- S + outer(sq, sq, "+")
  S[S < 0] <- 0
  K <- 1 - S / (2 * m)
  diag(K) <- 1
  K
}

