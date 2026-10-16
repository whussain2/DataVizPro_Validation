# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - ld_heatmap
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/genomicvisuals/linkage_disequilibrium/ld_server.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

ld_dprime_em <- function(M, max_iter = 100L, tol = 1e-7) {
  tryCatch({
    G <- round(M); G[!is.finite(G) | G < 0 | G > 2] <- NA
    m <- ncol(G)
    I0 <- (G == 0); I1 <- (G == 1); I2 <- (G == 2)
    I0[is.na(I0)] <- FALSE; I1[is.na(I1)] <- FALSE; I2[is.na(I2)] <- FALSE
    storage.mode(I0) <- "double"; storage.mode(I1) <- "double"; storage.mode(I2) <- "double"
    cnt <- function(A, B) crossprod(A, B)             # pairs with (a at i, b at j)
    n00 <- cnt(I0, I0); n01 <- cnt(I0, I1); n02 <- cnt(I0, I2)
    n10 <- cnt(I1, I0); n11 <- cnt(I1, I1); n12 <- cnt(I1, I2)
    n20 <- cnt(I2, I0); n21 <- cnt(I2, I1); n22 <- cnt(I2, I2)
    N <- n00 + n01 + n02 + n10 + n11 + n12 + n20 + n21 + n22
    # unambiguous haplotype counts (allele 1 = alternate)
    c00 <- 2 * n00 + n01 + n10; c01 <- 2 * n02 + n01 + n12
    c10 <- 2 * n20 + n10 + n21; c11 <- 2 * n22 + n21 + n12
    tot <- 2 * N; tot[tot == 0] <- 1
    f00 <- (c00 + 0.5 * n11) / tot; f01 <- (c01 + 0.5 * n11) / tot
    f10 <- (c10 + 0.5 * n11) / tot; f11 <- (c11 + 0.5 * n11) / tot
    for (it in seq_len(max_iter)) {
      den <- f00 * f11 + f01 * f10; den[den <= 0] <- 1
      cpl <- f00 * f11 / den                           # P(coupling | double het)
      cpl[f00 * f11 + f01 * f10 <= 0] <- 0.5
      g00 <- (c00 + cpl * n11) / tot; g01 <- (c01 + (1 - cpl) * n11) / tot
      g10 <- (c10 + (1 - cpl) * n11) / tot; g11 <- (c11 + cpl * n11) / tot
      delta <- max(abs(g00 - f00), abs(g01 - f01), abs(g10 - f10), abs(g11 - f11), na.rm = TRUE)
      f00 <- g00; f01 <- g01; f10 <- g10; f11 <- g11
      if (!is.finite(delta) || delta < tol) break
    }
    pA <- f10 + f11; pB <- f01 + f11                    # alternate-allele frequencies
    D  <- f11 - pA * pB
    Dmax_pos <- pmin(pA * (1 - pB), (1 - pA) * pB)
    Dmax_neg <- pmin(pA * pB, (1 - pA) * (1 - pB))
    Dmax <- ifelse(D >= 0, Dmax_pos, Dmax_neg)
    out <- ifelse(is.finite(Dmax) & Dmax > 0, abs(D) / Dmax, 0)
    out[!is.finite(out)] <- 0
    out <- pmin(pmax(out, 0), 1)
    out <- (out + t(out)) / 2
    diag(out) <- 1
    dimnames(out) <- list(colnames(M), colnames(M))
    out
  }, error = function(e) NULL)
}

ld_compute <- function(geno, measure = "r2") {
  M <- as.matrix(geno)
  storage.mode(M) <- "double"

  # Drop monomorphic markers - LD is undefined when a locus has no variation.
  keep <- apply(M, 2, function(x) {
    xx <- x[is.finite(x)]
    length(xx) >= 2 && stats::sd(xx) > 0
  })
  if (sum(keep) < 2) {
    stop("At least two polymorphic markers are required to compute LD. ",
         "After dropping monomorphic markers, fewer than two remain.")
  }
  M <- M[, keep, drop = FALSE]

  # Allelic correlation between marker columns (composite-LD estimator).
  R <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  R[!is.finite(R)] <- 0
  diag(R) <- 1

  if (identical(measure, "dprime")) {
    # HAP-005: D' from two-locus EM haplotype frequencies (Hill 1974 /
    # Excoffier & Slatkin 1995), as Haploview computes it, instead of the
    # composite-correlation approximation.
    LD <- ld_dprime_em(M)
    if (is.null(LD)) {
      LD <- matrix(0, ncol(M), ncol(M), dimnames = dimnames(R)); diag(LD) <- 1
    }
  } else {
    LD <- R^2
  }

  LD[!is.finite(LD)] <- 0
  LD <- pmin(pmax(LD, 0), 1)
  diag(LD) <- 1
  rownames(LD) <- colnames(LD) <- colnames(M)
  LD
}

