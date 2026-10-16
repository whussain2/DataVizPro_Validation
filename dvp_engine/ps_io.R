# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - ps_io
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/genomicvisuals/population_structures/ps_io.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

.ps_hwe_pvalues <- function(geno) {
  n_aa <- colSums(geno == 0, na.rm = TRUE)
  n_ab <- colSums(geno == 1, na.rm = TRUE)
  n_bb <- colSums(geno == 2, na.rm = TRUE)
  total <- n_aa + n_ab + n_bb
  freq_b <- (n_ab + 2 * n_bb) / (2 * total)
  freq_a <- 1 - freq_b
  exp_aa <- total * freq_a^2
  exp_ab <- total * 2 * freq_a * freq_b
  exp_bb <- total * freq_b^2
  chi <- (n_aa - exp_aa)^2 / exp_aa + (n_ab - exp_ab)^2 / exp_ab +
    (n_bb - exp_bb)^2 / exp_bb
  pvalue <- stats::pchisq(chi, df = 1L, lower.tail = FALSE)
  invalid <- total < 4 | freq_b <= 0 | freq_b >= 1 |
    exp_aa <= 0 | exp_ab <= 0 | exp_bb <= 0 | !is.finite(chi)
  pvalue[invalid] <- NA_real_
  pvalue
}

.ps_ld_prune <- function(geno, r2_threshold = 0.9, window = 50L) {
  imputed <- geno
  means <- colMeans(imputed, na.rm = TRUE); means[!is.finite(means)] <- 0
  na_cells <- which(!is.finite(imputed), arr.ind = TRUE)
  if (nrow(na_cells)) imputed[na_cells] <- means[na_cells[, 2L]]
  n_markers <- ncol(imputed)
  window <- max(2L, min(as.integer(window), n_markers))
  centered <- imputed - rep(colMeans(imputed), each = nrow(imputed))
  norms <- sqrt(colSums(centered^2)); norms[!is.finite(norms) | norms == 0] <- 1
  standardized <- centered / rep(norms, each = nrow(imputed))
  threshold <- sqrt(r2_threshold)
  keep <- logical(n_markers)
  retained_recent <- integer(0)
  for (index in seq_len(n_markers)) {
    drop_marker <- FALSE
    if (length(retained_recent)) {
      correlations <- crossprod(standardized[, index],
        standardized[, retained_recent, drop = FALSE])
      if (any(abs(correlations) > threshold, na.rm = TRUE)) drop_marker <- TRUE
    }
    if (!drop_marker) {
      keep[index] <- TRUE
      retained_recent <- c(retained_recent, index)
      if (length(retained_recent) > window) retained_recent <- retained_recent[-1L]
    }
  }
  keep
}

ps_prepare_matrix <- function(geno, ids = rownames(geno), groups = NULL, metadata = NULL,
                              marker_map = NULL, source = NULL, min_sample_call = 0.8,
                              min_marker_call = 0.8, min_maf = 0.01,
                              max_maf = 0.5, max_markers = Inf, seed = 42L,
                              drop_duplicates = FALSE, hwe_p = 0,
                              ld_prune = FALSE, ld_r2 = 0.9, ld_window = 50L) {
  geno <- as.matrix(geno)
  storage.mode(geno) <- "double"
  if (is.null(ids)) ids <- paste0("Sample_", seq_len(nrow(geno)))
  ids <- make.unique(as.character(ids))
  rownames(geno) <- ids
  geno[!is.finite(geno)] <- NA_real_
  if (any(geno < 0 | geno > 2, na.rm = TRUE)) stop("Prepared genotypes must use 0/1/2 allele counts.")
  if (any(abs(geno - round(geno)) > sqrt(.Machine$double.eps), na.rm = TRUE)) {
    stop("Prepared genotypes must contain discrete 0/1/2 allele counts; fractional values cannot be analyzed consistently by sNMF and the genetic tree.")
  }
  geno <- round(geno)

  sample_call <- rowMeans(!is.na(geno))
  keep_samples <- is.finite(sample_call) & sample_call >= min_sample_call
  geno <- geno[keep_samples, , drop = FALSE]
  ids <- ids[keep_samples]
  if (!is.null(groups)) groups <- groups[keep_samples]
  if (!is.null(metadata)) metadata <- metadata[keep_samples, , drop = FALSE]

  # Vectorized marker QC (colSums/colMeans only) — an apply() var scan over
  # hundreds of thousands of markers is the main preparation bottleneck.
  n_obs <- colSums(!is.na(geno))
  marker_call <- n_obs / nrow(geno)
  col_mean <- colSums(geno, na.rm = TRUE) / pmax(n_obs, 1)
  frequency <- col_mean / 2
  maf <- pmin(frequency, 1 - frequency)
  variance <- (colSums(geno * geno, na.rm = TRUE) - n_obs * col_mean^2) / pmax(n_obs - 1, 1)
  variance[n_obs < 2] <- NA_real_
  keep_markers <- is.finite(marker_call) & marker_call >= min_marker_call &
    is.finite(maf) & maf >= min_maf & maf <= max_maf &
    is.finite(variance) & variance > 0
  hwe_pvalue <- rep(NA_real_, ncol(geno))
  if (is.finite(hwe_p) && hwe_p > 0) {
    hwe_pvalue <- .ps_hwe_pvalues(geno)
    keep_markers <- keep_markers & (!is.finite(hwe_pvalue) | hwe_pvalue >= hwe_p)
  }
  marker_type <- if (is.null(source$marker_type)) "" else source$marker_type
  microsatellite <- grepl("microsatellite", marker_type, fixed = TRUE)
  qc_notes <- character(0)
  if (microsatellite && !is.null(marker_map) && "locus" %in% names(marker_map)) {
    # POP-004: microsatellite filters work per LOCUS (heterozygosity and the
    # tree distance need each locus's allele set): a locus is kept when one
    # of its alleles passes the call-rate / variance filters AND its locus
    # MAF (1 - frequency of its most common allele) reaches Minimum MAF.
    locus_chr <- as.character(marker_map$locus)
    allele_ok <- is.finite(marker_call) & marker_call >= min_marker_call &
      is.finite(variance) & variance > 0
    if (is.finite(hwe_p) && hwe_p > 0) allele_ok <- allele_ok & (!is.finite(hwe_pvalue) | hwe_pvalue >= hwe_p)
    major <- tapply(ifelse(is.finite(frequency), frequency, 0), locus_chr, max)
    locus_maf <- 1 - major
    eligible_loci <- intersect(unique(locus_chr[allele_ok]),
                               names(locus_maf)[locus_maf >= min_maf])
    # constant alleles (never or always observed) carry no information and
    # cannot be scaled for PCA, so they are dropped even inside kept loci
    keep_markers <- locus_chr %in% eligible_loci & is.finite(variance) & variance > 0
    qc_notes <- c(qc_notes, sprintf(
      "Microsatellites: QC is applied per locus (%d of %d loci kept, with their informative alleles); Minimum MAF uses 1 - the most common allele's frequency.%s",
      length(eligible_loci), length(unique(locus_chr)),
      if (isTRUE(ld_prune)) " LD pruning does not apply to multi-allelic loci and was skipped." else ""))
  }
  qc <- data.frame(marker = colnames(geno), call_rate = marker_call, maf = maf,
    variance = variance, hwe_p = hwe_pvalue, retained = keep_markers,
    stringsAsFactors = FALSE)
  geno <- geno[, keep_markers, drop = FALSE]
  if (!is.null(marker_map) && "marker" %in% names(marker_map)) {
    marker_map <- marker_map[match(colnames(geno), marker_map$marker), , drop = FALSE]
  }
  if (nrow(geno) < 4 || ncol(geno) < 2) stop("Too few samples or polymorphic markers remain after QC.")
  # Cap the marker count BEFORE the O(markers x window) duplicate and LD passes
  # so they only ever operate on the budgeted subset — this keeps preparation
  # tractable for panels with hundreds of thousands or millions of markers.
  if (is.finite(max_markers) && ncol(geno) > max_markers) {
    set.seed(seed)
    if (microsatellite && !is.null(marker_map) && "locus" %in% names(marker_map)) {
      loci <- sample(unique(as.character(marker_map$locus)))
      locus_sizes <- table(as.character(marker_map$locus))[loci]
      selected_loci <- loci[cumsum(locus_sizes) <= max_markers]
      if (!length(selected_loci)) selected_loci <- loci[[1]]
      selected <- which(as.character(marker_map$locus) %in% selected_loci)
    } else {
      selected <- sort(sample(seq_len(ncol(geno)), as.integer(max_markers)))
    }
    geno <- geno[, selected, drop = FALSE]
    if (!is.null(marker_map)) marker_map <- marker_map[selected, , drop = FALSE]
  }
  if (isTRUE(drop_duplicates) && ncol(geno) > 1L) {
    duplicated_markers <- duplicated(t(geno))
    geno <- geno[, !duplicated_markers, drop = FALSE]
    if (!is.null(marker_map)) marker_map <- marker_map[!duplicated_markers, , drop = FALSE]
  }
  if (isTRUE(ld_prune) && ncol(geno) > 1L && !microsatellite) {
    keep_ld <- .ps_ld_prune(geno, r2_threshold = ld_r2, window = as.integer(ld_window))
    geno <- geno[, keep_ld, drop = FALSE]
    if (!is.null(marker_map)) marker_map <- marker_map[keep_ld, , drop = FALSE]
  }
  if (ncol(geno) < 2) stop("Too few polymorphic markers remain after de-duplication and LD pruning.")
  list(geno = geno, ids = ids, groups = groups, metadata = metadata,
    marker_map = marker_map, qc = qc, source = source,
    n_input_samples = length(keep_samples), n_input_markers = length(marker_call),
    sample_call_rate = sample_call, pre_qc = TRUE, qc_notes = qc_notes)
}

