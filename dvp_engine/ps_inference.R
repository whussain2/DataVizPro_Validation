# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - ps_inference
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/genomicvisuals/population_structures/ps_inference.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

.ps_impute_genotypes <- function(geno) {
  result <- geno
  means <- colMeans(result, na.rm = TRUE)
  means[!is.finite(means)] <- 0
  missing <- which(!is.finite(result), arr.ind = TRUE)
  if (nrow(missing)) result[missing] <- means[missing[, 2]]
  result
}

ps_run_pca <- function(prepared, components = 10L, scalable = FALSE) {
  matrix <- .ps_impute_genotypes(prepared$geno)
  components <- max(2L, min(as.integer(components), min(dim(matrix)) - 1L))
  use_irlba <- isTRUE(scalable) && requireNamespace("irlba", quietly = TRUE) &&
    min(dim(matrix)) > components + 1L
  fit <- if (use_irlba) {
    irlba::prcomp_irlba(matrix, n = components, center = TRUE, scale. = TRUE)
  } else {
    stats::prcomp(matrix, center = TRUE, scale. = TRUE, rank. = components)
  }
  # POP-002: share of the TOTAL variance. After scaling every marker has
  # variance 1, so the total is the number of markers - not the sum of the
  # few components a truncated (irlba) PCA returns
  variance <- fit$sdev^2 / ncol(matrix)
  scores <- as.data.frame(fit$x, check.names = FALSE)
  scores$sample_id <- prepared$ids
  if (!is.null(prepared$groups)) scores$known_group <- prepared$groups
  list(fit = fit, scores = scores, variance = variance,
       engine = if (use_irlba) "irlba truncated PCA" else "stats::prcomp")
}

.ps_admix_clamp <- function(value, epsilon = 1e-8) {
  pmin(pmax(value, epsilon), 1 - epsilon)
}

.ps_admix_block_size <- function(n_individuals, target_cells = 2e6) {
  max(1L, floor(target_cells / max(1L, n_individuals)))
}

.ps_admix_initialize <- function(geno, K, seed) {
  set.seed(as.integer(seed))
  n <- nrow(geno)
  m <- ncol(geno)
  q <- matrix(stats::rexp(n * K), nrow = n, ncol = K)
  q <- q / rowSums(q)
  frequency <- colMeans(geno, na.rm = TRUE) / 2
  frequency[!is.finite(frequency)] <- 0.5
  frequency <- .ps_admix_clamp(frequency, 1e-4)
  p <- matrix(rep(stats::qlogis(frequency), each = K), nrow = K)
  if (K > 1L) p <- p + matrix(stats::rnorm(K * m, sd = 0.25), nrow = K)
  p <- stats::plogis(p)
  list(q = q, p = .ps_admix_clamp(p))
}

.ps_admix_loglik <- function(geno, q, p) {
  block <- .ps_admix_block_size(nrow(geno))
  total <- 0
  for (start in seq.int(1L, ncol(geno), by = block)) {
    columns <- start:min(start + block - 1L, ncol(geno))
    observed <- is.finite(geno[, columns, drop = FALSE])
    dosage <- geno[, columns, drop = FALSE]
    dosage[!observed] <- 0
    frequency <- .ps_admix_clamp(q %*% p[, columns, drop = FALSE])
    total <- total + sum(dosage * log(frequency) +
      (2 * observed - dosage) * log1p(-frequency))
  }
  total
}

.ps_run_em_fit <- function(geno, K, iterations, tolerance, seed) {
  state <- .ps_admix_initialize(geno, K, seed)
  q <- state$q
  p <- state$p
  block <- .ps_admix_block_size(nrow(geno))
  observed_count <- max(1L, sum(is.finite(geno)))
  previous <- -Inf
  converged <- FALSE
  iteration <- 0L

  for (iteration in seq_len(max(1L, as.integer(iterations)))) {
    q_numerator <- matrix(0, nrow(geno), K)
    p_numerator <- matrix(0, K, ncol(geno))
    p_denominator <- matrix(0, K, ncol(geno))
    log_likelihood <- 0

    for (start in seq.int(1L, ncol(geno), by = block)) {
      columns <- start:min(start + block - 1L, ncol(geno))
      values <- geno[, columns, drop = FALSE]
      observed <- is.finite(values)
      values[!observed] <- 0
      alternate <- 2 * observed - values
      p_block <- p[, columns, drop = FALSE]
      frequency <- .ps_admix_clamp(q %*% p_block)
      log_likelihood <- log_likelihood + sum(values * log(frequency) +
        alternate * log1p(-frequency))

      for (cluster in seq_len(K)) {
        reference_weight <- outer(q[, cluster], p_block[cluster, ]) / frequency
        alternate_weight <- outer(q[, cluster], 1 - p_block[cluster, ]) /
          (1 - frequency)
        reference_count <- values * reference_weight
        alternate_count <- alternate * alternate_weight
        q_numerator[, cluster] <- q_numerator[, cluster] +
          rowSums(reference_count) + rowSums(alternate_count)
        p_numerator[cluster, columns] <- colSums(reference_count)
        p_denominator[cluster, columns] <- colSums(reference_count) +
          colSums(alternate_count)
      }
    }

    q_total <- rowSums(q_numerator)
    q_total[q_total <= 0 | !is.finite(q_total)] <- 1
    q <- q_numerator / q_total
    updated_p <- p_numerator / pmax(p_denominator, 1e-12)
    invalid <- !is.finite(updated_p)
    updated_p[invalid] <- p[invalid]
    p <- .ps_admix_clamp(updated_p)

    if (is.finite(previous) &&
        abs(log_likelihood - previous) / observed_count < tolerance) {
      converged <- TRUE
      break
    }
    previous <- log_likelihood
  }

  list(q = q, p = p, log_likelihood = .ps_admix_loglik(geno, q, p),
       iterations = iteration, converged = converged)
}

.ps_run_vb_fit <- function(geno, K, iterations, tolerance, seed,
                           alpha = 1, beta = 1) {
  state <- .ps_admix_initialize(geno, K, seed)
  q <- state$q
  p <- state$p
  block <- .ps_admix_block_size(nrow(geno))
  observed_count <- max(1L, sum(is.finite(geno)))
  individual_observed <- rowSums(is.finite(geno))
  marker_observed <- colSums(is.finite(geno))
  q_shape <- alpha + q * (2 * individual_observed)
  p_reference <- beta + p * matrix(2 * marker_observed, K, ncol(geno), byrow = TRUE)
  p_alternate <- beta + (1 - p) *
    matrix(2 * marker_observed, K, ncol(geno), byrow = TRUE)
  previous <- -Inf
  converged <- FALSE
  iteration <- 0L

  for (iteration in seq_len(max(1L, as.integer(iterations)))) {
    expected_q <- exp(digamma(pmax(q_shape, 1e-12)) -
      digamma(pmax(rowSums(q_shape), 1e-12)))
    p_total <- p_reference + p_alternate
    expected_reference <- exp(digamma(pmax(p_reference, 1e-12)) -
      digamma(pmax(p_total, 1e-12)))
    expected_alternate <- exp(digamma(pmax(p_alternate, 1e-12)) -
      digamma(pmax(p_total, 1e-12)))
    q_numerator <- matrix(0, nrow(geno), K)
    reference_numerator <- matrix(0, K, ncol(geno))
    alternate_numerator <- matrix(0, K, ncol(geno))
    log_likelihood <- 0

    for (start in seq.int(1L, ncol(geno), by = block)) {
      columns <- start:min(start + block - 1L, ncol(geno))
      values <- geno[, columns, drop = FALSE]
      observed <- is.finite(values)
      values[!observed] <- 0
      alternate <- 2 * observed - values
      frequency <- .ps_admix_clamp(q %*% p[, columns, drop = FALSE])
      log_likelihood <- log_likelihood + sum(values * log(frequency) +
        alternate * log1p(-frequency))
      expected_reference_block <- expected_reference[, columns, drop = FALSE]
      expected_alternate_block <- expected_alternate[, columns, drop = FALSE]
      reference_total <- pmax(expected_q %*% expected_reference_block, 1e-300)
      alternate_total <- pmax(expected_q %*% expected_alternate_block, 1e-300)

      for (cluster in seq_len(K)) {
        reference_weight <- outer(expected_q[, cluster],
          expected_reference_block[cluster, ]) / reference_total
        alternate_weight <- outer(expected_q[, cluster],
          expected_alternate_block[cluster, ]) / alternate_total
        reference_count <- values * reference_weight
        alternate_count <- alternate * alternate_weight
        q_numerator[, cluster] <- q_numerator[, cluster] +
          rowSums(reference_count) + rowSums(alternate_count)
        reference_numerator[cluster, columns] <- colSums(reference_count)
        alternate_numerator[cluster, columns] <- colSums(alternate_count)
      }
    }

    q_shape <- alpha + q_numerator
    p_reference <- beta + reference_numerator
    p_alternate <- beta + alternate_numerator
    q <- q_shape / rowSums(q_shape)
    p <- .ps_admix_clamp(p_reference / (p_reference + p_alternate))

    if (is.finite(previous) &&
        abs(log_likelihood - previous) / observed_count < tolerance) {
      converged <- TRUE
      break
    }
    previous <- log_likelihood
  }

  list(q = q, p = p, log_likelihood = .ps_admix_loglik(geno, q, p),
       iterations = iteration, converged = converged)
}

ps_run_parametric_admixture <- function(prepared, engine = c("em", "vb"),
                                        k_min = 1L, k_max = 8L,
                                        repetitions = 5L,
                                        entropy_fraction = 0.05,
                                        iterations = 200L, seed = 42L,
                                        tolerance = 0.0001,
                                        vb_alpha = 1, vb_beta = 1) {
  engine <- match.arg(engine)
  if (!is.finite(iterations) || iterations < 1L)
    stop("Iterations must be a positive integer.")
  if (!is.finite(repetitions) || repetitions < 1L)
    stop("Runs per K must be a positive integer.")
  if (!is.finite(entropy_fraction) || entropy_fraction <= 0 || entropy_fraction >= 1)
    stop("Masked genotype fraction must be greater than 0 and less than 1.")
  if (!is.finite(tolerance) || tolerance <= 0)
    stop("Convergence tolerance must be positive.")
  if (!is.finite(seed)) stop("Random seed must be finite.")
  if (identical(engine, "vb") &&
      (!is.finite(vb_alpha) || vb_alpha <= 0 || !is.finite(vb_beta) || vb_beta <= 0))
    stop("Variational Bayes alpha and beta priors must be positive.")
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue is required to align ancestry labels across runs.")
  }
  geno <- as.matrix(prepared$geno)
  storage.mode(geno) <- "double"
  if (!all(geno[is.finite(geno)] %in% 0:2)) {
    stop("EM and variational Bayes require discrete 0/1/2 genotype dosages.")
  }
  if (ncol(geno) > 20000L) {
    stop(paste0(
      "EM and variational Bayes are limited to 20,000 retained markers to keep ",
      "interactive runs reliable. Reduce Maximum markers or use LEA sNMF for genome-scale data."))
  }
  k_values <- seq.int(max(1L, as.integer(k_min)), max(2L, as.integer(k_max)))
  if (max(k_values) > nrow(geno)) {
    stop("The maximum K cannot exceed the number of retained individuals.")
  }
  repetitions <- max(1L, as.integer(repetitions))
  observed_indices <- which(is.finite(geno))
  if (length(observed_indices) < 100L) {
    stop("At least 100 observed genotype calls are required for masked validation.")
  }
  masked_count <- max(10L, floor(length(observed_indices) * entropy_fraction))
  masked_count <- min(masked_count, length(observed_indices) - 1L)
  masks <- lapply(seq_len(repetitions), function(run) {
    set.seed(as.integer(seed) + run * 100003L)
    sample(observed_indices, masked_count)
  })

  q_by_k <- list()
  run_table <- list()
  for (k in k_values) {
    fits <- lapply(seq_len(repetitions), function(run) {
      training <- geno
      held <- masks[[run]]
      truth <- training[held]
      training[held] <- NA_real_
      fit_seed <- as.integer(seed) + k * 1009L + run * 100003L
      fit <- if (identical(engine, "em")) {
        .ps_run_em_fit(training, k, iterations, tolerance, fit_seed)
      } else {
        .ps_run_vb_fit(training, k, iterations, tolerance, fit_seed,
          alpha = vb_alpha, beta = vb_beta)
      }
      location <- arrayInd(held, dim(geno))
      predicted_frequency <- rowSums(
        fit$q[location[, 1], , drop = FALSE] *
          t(fit$p[, location[, 2], drop = FALSE]))
      predicted_frequency <- .ps_admix_clamp(predicted_frequency)
      fit$cross_entropy <- -mean(truth * log(predicted_frequency) +
        (2 - truth) * log1p(-predicted_frequency)) / 2
      # Evaluate every replicate on the same complete observed-call set so
      # likelihood differences across K and runs are directly comparable.
      fit$log_likelihood <- .ps_admix_loglik(geno, fit$q, fit$p)
      rownames(fit$q) <- prepared$ids
      colnames(fit$q) <- paste0("Cluster ", seq_len(k))
      fit
    })
    entropy <- vapply(fits, `[[`, numeric(1), "cross_entropy")
    reference_run <- which.min(entropy)
    reference <- fits[[reference_run]]$q
    aligned <- lapply(seq_along(fits), function(run) {
      if (run == reference_run) list(q = reference, similarity = 1) else
        .ps_align_q(reference, fits[[run]]$q)
    })
    consensus <- Reduce("+", lapply(aligned, `[[`, "q")) / length(aligned)
    consensus <- consensus / rowSums(consensus)
    stability <- vapply(aligned, `[[`, numeric(1), "similarity")
    q_by_k[[as.character(k)]] <- list(
      consensus = consensus, best = reference,
      runs = lapply(aligned, `[[`, "q"), entropy = entropy,
      best_run = reference_run, stability = stability)
    run_table[[as.character(k)]] <- data.frame(
      K = k, run = seq_len(repetitions), cross_entropy = entropy,
      alignment_similarity = stability,
      log_likelihood = vapply(fits, `[[`, numeric(1), "log_likelihood"),
      iterations = vapply(fits, `[[`, integer(1), "iterations"),
      converged = vapply(fits, `[[`, logical(1), "converged"),
      stringsAsFactors = FALSE)
  }
  runs <- do.call(rbind, run_table)
  summary <- do.call(rbind, lapply(split(runs, runs$K), function(frame) {
    data.frame(K = frame$K[1], mean_cross_entropy = mean(frame$cross_entropy),
      sd_cross_entropy = stats::sd(frame$cross_entropy),
      best_cross_entropy = min(frame$cross_entropy),
      mean_stability = mean(frame$alignment_similarity, na.rm = TRUE),
      mean_log_likelihood = mean(frame$log_likelihood),
      sd_log_likelihood = stats::sd(frame$log_likelihood))
  }))
  rownames(summary) <- NULL
  best_k <- summary$K[which.min(summary$mean_cross_entropy)]
  k_selection <- .ps_k_selection(summary)
  delta_k <- k_selection$delta_k
  evanno_k <- if (any(is.finite(delta_k)))
    k_selection$K[which.max(replace(delta_k, !is.finite(delta_k), -Inf))] else NA_integer_
  engine_label <- if (identical(engine, "em"))
    "Native EM maximum-likelihood admixture" else
    "Native mean-field variational Bayes admixture"
  nonconverged_runs <- sum(!runs$converged)
  k_rationale <- .ps_k_rationale(best_k, evanno_k, engine)
  if (nonconverged_runs > 0L) k_rationale <- paste(k_rationale, sprintf(
    "%d of %d runs reached the iteration limit; increase Iterations or review those diagnostics before interpretation.",
    nonconverged_runs, nrow(runs)))
  list(q_by_k = q_by_k, runs = runs, summary = summary, suggested_k = best_k,
       k_selection = k_selection, suggested_k_entropy = best_k,
       suggested_k_evanno = evanno_k,
       k_rationale = k_rationale,
       settings = list(k_min = min(k_values), k_max = max(k_values),
         repetitions = repetitions, entropy_fraction = entropy_fraction,
         iterations = iterations, seed = seed, tolerance = tolerance,
         vb_alpha = vb_alpha, vb_beta = vb_beta),
      nonconverged_runs = nonconverged_runs,
       engine = engine_label, engine_id = engine)
}

.ps_write_lea_geno <- function(geno, path) {
  encoded <- t(geno)
  encoded[!is.finite(encoded)] <- 9
  encoded <- round(encoded)
  writeLines(apply(encoded, 1, paste0, collapse = ""), path)
  invisible(path)
}

.ps_align_q <- function(reference, candidate) {
  if (ncol(reference) != ncol(candidate)) stop("Q matrices must have the same K.")
  denominator <- outer(sqrt(colSums(reference^2)), sqrt(colSums(candidate^2)))
  similarity <- crossprod(reference, candidate) / denominator
  similarity[!is.finite(similarity)] <- 0
  assignment <- clue::solve_LSAP(similarity, maximum = TRUE)
  aligned <- candidate[, as.integer(assignment), drop = FALSE]
  colnames(aligned) <- colnames(reference)
  aligned_denominator <- sqrt(colSums(reference^2)) * sqrt(colSums(aligned^2))
  aligned_similarity <- colSums(reference * aligned) / aligned_denominator
  list(q = aligned, similarity = mean(aligned_similarity, na.rm = TRUE))
}

.ps_k_selection <- function(summary) {
  sel <- summary[order(summary$K), , drop = FALSE]
  rownames(sel) <- NULL
  n <- nrow(sel)
  has_likelihood <- all(c("mean_log_likelihood", "sd_log_likelihood") %in% names(sel))
  lk <- if (has_likelihood) sel$mean_log_likelihood else rep(NA_real_, n)
  sdk <- if (has_likelihood) sel$sd_log_likelihood else rep(NA_real_, n)
  l_prime <- rep(NA_real_, n)
  l_double_prime <- rep(NA_real_, n)
  delta_k <- rep(NA_real_, n)
  if (n >= 2L) for (i in 2:n) l_prime[i] <- lk[i] - lk[i - 1L]
  if (n >= 3L) for (i in 2:(n - 1L)) {
    l_double_prime[i] <- abs(lk[i + 1L] - 2 * lk[i] + lk[i - 1L])
    if (is.finite(sdk[i]) && sdk[i] > 0) delta_k[i] <- l_double_prime[i] / sdk[i]
  }
  sel$L_k <- lk
  sel$L_prime <- l_prime
  sel$L_double_prime <- l_double_prime
  sel$delta_k <- delta_k
  sel
}

.ps_k_rationale <- function(entropy_k, evanno_k, engine = "snmf") {
  parts <- if (identical(engine, "snmf")) sprintf(
    "Cross-entropy is minimized at K = %d (LEA sNMF criterion; Frichot & Fran\u00e7ois 2015).",
    entropy_k) else sprintf(
    "Masked-genotype cross-entropy is minimized at K = %d for the selected inference engine.",
    entropy_k)
  if (!is.na(evanno_k)) {
    parts <- c(parts, sprintf(
      "The Evanno-style \u0394K rate-of-change peaks at K = %d (after Evanno et al. 2005), which may mark an upper hierarchical level of structure.",
      evanno_k))
    parts <- c(parts, if (identical(as.integer(evanno_k), as.integer(entropy_k)))
      "Both criteria agree, giving strong support for this K." else
      "When the criteria disagree, favour the smaller K unless a larger K is biologically meaningful; \u0394K cannot evaluate K = 1 (no structure).")
  } else {
    parts <- c(parts, if (identical(engine, "snmf"))
      "Evanno-style delta-K is not reported for sNMF because cross-entropy is predictive error, not replicate log-likelihood." else
      "The Evanno-style delta-K statistic needs at least three K values with replicate log-likelihood variance and is not available for this run.")
  }
  paste(parts, collapse = " ")
}

ps_run_snmf <- function(prepared, k_min = 1L, k_max = 8L, repetitions = 5L,
                        entropy_fraction = 0.05, alpha = 100, iterations = 200L,
                        cpu = 1L, seed = 42L, tolerance = 0.0001) {
  if (!requireNamespace("LEA", quietly = TRUE)) {
    stop("LEA is required for sNMF ancestry inference.")
  }
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue is required to align ancestry labels across runs.")
  }
  k_values <- seq.int(max(1L, as.integer(k_min)), max(2L, as.integer(k_max)))
  repetitions <- max(1L, as.integer(repetitions))
  path <- tempfile("datavizpro_structure_", fileext = ".geno")
  .ps_write_lea_geno(prepared$geno, path)
  base <- tools::file_path_sans_ext(path)
  on.exit(unlink(c(path, paste0(base, ".snmfProject"), paste0(base, ".snmf")),
                 recursive = TRUE, force = TRUE), add = TRUE)
  set.seed(seed)
  output <- utils::capture.output({
    project <- LEA::snmf(path, K = k_values, repetitions = repetitions,
      CPU = max(1L, as.integer(cpu)), alpha = alpha, tolerance = tolerance,
      entropy = TRUE, percentage = entropy_fraction,
      iterations = as.integer(iterations), project = "new", seed = as.integer(seed))
  })

  q_by_k <- list()
  run_table <- list()
  for (k in k_values) {
    entropy <- as.numeric(LEA::cross.entropy(project, K = k))
    q_runs <- lapply(seq_len(repetitions), function(run) {
      q <- LEA::Q(project, K = k, run = run)
      rownames(q) <- prepared$ids
      colnames(q) <- paste0("Cluster ", seq_len(k))
      q
    })
    reference_run <- which.min(entropy)
    reference <- q_runs[[reference_run]]
    aligned <- lapply(seq_along(q_runs), function(run) {
      if (run == reference_run) list(q = reference, similarity = 1) else
        .ps_align_q(reference, q_runs[[run]])
    })
    consensus <- Reduce("+", lapply(aligned, `[[`, "q")) / length(aligned)
    consensus <- consensus / rowSums(consensus)
    stability <- vapply(aligned, `[[`, numeric(1), "similarity")
    q_by_k[[as.character(k)]] <- list(
      consensus = consensus, best = reference, runs = lapply(aligned, `[[`, "q"),
      entropy = entropy, best_run = reference_run, stability = stability
    )
    run_table[[as.character(k)]] <- data.frame(
      K = k, run = seq_len(repetitions), cross_entropy = entropy,
      alignment_similarity = stability, stringsAsFactors = FALSE)
  }
  runs <- do.call(rbind, run_table)
  summary <- do.call(rbind, lapply(split(runs, runs$K), function(frame) {
    data.frame(K = frame$K[1], mean_cross_entropy = mean(frame$cross_entropy),
      sd_cross_entropy = stats::sd(frame$cross_entropy),
      best_cross_entropy = min(frame$cross_entropy),
      mean_stability = mean(frame$alignment_similarity, na.rm = TRUE))
  }))
  rownames(summary) <- NULL
  best_k <- summary$K[which.min(summary$mean_cross_entropy)]
  k_selection <- .ps_k_selection(summary)
  entropy_k <- best_k
  delta_k <- k_selection$delta_k
  evanno_k <- if (any(is.finite(delta_k)))
    k_selection$K[which.max(replace(delta_k, !is.finite(delta_k), -Inf))] else NA_integer_
  k_rationale <- .ps_k_rationale(entropy_k, evanno_k, "snmf")
  list(q_by_k = q_by_k, runs = runs, summary = summary, suggested_k = best_k,
       k_selection = k_selection, suggested_k_entropy = entropy_k,
       suggested_k_evanno = evanno_k, k_rationale = k_rationale,
       settings = list(k_min = min(k_values), k_max = max(k_values),
         repetitions = repetitions, entropy_fraction = entropy_fraction,
         alpha = alpha, iterations = iterations, cpu = cpu, seed = seed,
         tolerance = tolerance),
      log_tail = utils::tail(output, 20), engine = "LEA::snmf", engine_id = "snmf")
}

