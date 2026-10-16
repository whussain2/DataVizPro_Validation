# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - ld_decay
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/genomicvisuals/ld_decay/ldd_analysis.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

ldd_default_settings <- function() {
  list(
    view_mode = "decay",
    group_by = "none",              # none | chromosome | population (set in Colours panel)
    facet_enabled = FALSE,          # separate panels per group (Multi-Panel panel)
    density = FALSE,                # density heat-scatter override (Appearance panel)
    model = "hw",                   # hw (Hill-Weir) | exp (Sved/exponential) | none
    threshold_r2 = 0.1,             # practical "extent of LD" cutoff (fixed value)
    threshold_source = "fixed",     # fixed | background95 (95th pct of unlinked r^2)
    max_dist_display = 0,           # 0 = auto (use full window)

    # LD calculation controls (Data Preparation)
    ld_measure = "r2",              # r2 (composite squared correlation) | dprime
    window_mode = "distance",       # distance (bp) | snp_count | both
    max_snp_window = 100L,          # SNPs on each side when window_mode uses counts
    min_r2_keep = 0,                # drop pairs with r^2 below this (PLINK --ld-window-r2)
    max_het = 1,                    # max heterozygosity ratio per marker (1 = off)
    pair_cap = 2e5,                 # random-subsample pair cap per chromosome
    pair_seed = 1L,                 # seed for reproducible pair subsampling

    # Distance binning controls
    bin_mode = "auto",              # auto | fixed_width | n_bins | log
    bin_width_kb = 0,               # bin width (kb) when fixed_width
    n_bins = 200L,                  # target bin count when n_bins / log

    # Titles
    title = "Linkage Disequilibrium Decay",
    subtitle = "",
    title_size = 16, title_color = "#0f172a", title_family = "sans",
    title_face = "bold", title_hjust = 0.5,
    subtitle_size = 12, subtitle_color = "#475569", subtitle_family = "sans",
    subtitle_face = "italic",

    # Axis titles (overrides per view via *_label_overrides)
    x_label_overrides = list(), y_label_overrides = list(),
    x_title_size = 12, x_title_color = "#0f172a", x_title_family = "sans",
    x_title_face = "bold", x_title_hide = FALSE,
    y_title_size = 12, y_title_color = "#0f172a", y_title_family = "sans",
    y_title_face = "bold", y_title_hide = FALSE,

    # Axis tick text (X and Y styled independently)
    axis_text_size = 10, axis_text_color = "#334155", axis_text_family = "sans",
    axis_text_face = "plain", x_angle = 0,
    x_tick_size = 10, x_tick_color = "#334155", x_tick_family = "sans", x_tick_face = "plain",
    y_tick_size = 10, y_tick_color = "#334155", y_tick_family = "sans", y_tick_face = "plain",
    x_tick_length = 0.15, x_tick_line_color = "#334155", x_tick_linetype = "solid", x_tick_linewidth = 0.5,
    y_tick_length = 0.15, y_tick_line_color = "#334155", y_tick_linetype = "solid", y_tick_linewidth = 0.5,
    x_scale = "linear",             # linear | log10 (distance axis)
    dist_unit = "auto",             # auto | bp | kb | Mb (axis label formatting)

    # Appearance
    show_points = TRUE, point_color = "#94a3b8", point_size = 0.7, point_alpha = 0.25,
    show_binned = TRUE, binned_color = "#2C7FB8", binned_size = 1.4, binned_alpha = 0.8,
    show_fit = TRUE, fit_color = "#D7301F", fit_width = 1.3, fit_linetype = "solid",
    show_ci = FALSE, ci_alpha = 0.15, ci_boot = 200L,
    show_half_decay = TRUE, show_threshold = TRUE, guide_color = "#64748b",
    show_metrics_box = TRUE, metrics_box_corner = "topright",
    metrics_box_size = 3.2, metrics_box_text_color = "#1e293b",
    metrics_box_fill = "#ffffff", metrics_box_fill_alpha = 0.7, metrics_box_face = "plain",
    metrics_box_show_half_decay = TRUE, metrics_box_half_decay_label = "Half-decay",
    metrics_box_show_extent = TRUE, metrics_box_extent_label = "LD extent",
    metrics_box_show_intercept = TRUE, metrics_box_intercept_label = "Intercept~italic(r)^2",

    # Colours (grouped views)
    color_mode = "palette", palette = "Set2", manual_colors = list(),
    reverse_palette = FALSE,

    # Density view
    density_bins = 60, density_bins_x = 80L, density_bins_y = 50L,
    density_measure = "count", density_palette = "viridis", density_reverse = FALSE,
    density_trans = "log10", density_clip_quantile = 1,
    density_show_contours = FALSE, density_contour_bins = 5L,
    density_contour_color = "#ffffff", density_contour_alpha = 0.65,
    density_profile = "none", density_profile_color = "#ffffff",
    density_profile_width = 1.1, density_legend_show = TRUE,

    # Metrics (bar) view
    metric_stat = "half_decay",     # half_decay | threshold
    bar_sort = "value",             # value | name

    # Theme
    base_size = 13, theme = "classic", background = "#ffffff",
    panel_background = "#ffffff", grid_show = FALSE, grid_color = "#e5e7eb",
    grid_minor_show = FALSE,

    # Legend
    legend_show = TRUE, legend_position = "right", legend_title = "Group",
    legend_text_size = 11, legend_text_color = "#334155",
    legend_text_family = "sans", legend_text_face = "plain",
    legend_title_size = 12, legend_title_color = "#0f172a",
    legend_title_family = "sans", legend_title_face = "bold", legend_title_hjust = 0,
    legend_ncol = 1L, legend_reverse = FALSE, legend_direction = "vertical",
    legend_key_width = 6, legend_key_height = 6, legend_spacing = 4,
    legend_inside_x = 0.95, legend_inside_y = 0.95,
    legend_bg_show = FALSE, legend_bg_color = "#ffffff", legend_bg_alpha = 0.8,
    legend_border_color = "#cccccc", legend_border_width = 0,
    legend_labels = list(),

    # Facet
    facet_ncol = 3L, facet_scales = "fixed", facet_strip_size = 11,
    facet_strip_color = "#0f172a", facet_strip_bg = "#e2e8f0",
    facet_strip_face = "bold", facet_strip_border = "#cbd5e1", facet_strip_border_width = 0,

    # Annotations
    annotations_by_view = list()
  )
}

ldd_maf <- function(X) {
  p <- colMeans(X, na.rm = TRUE) / 2
  pmin(p, 1 - p)
}

ldd_filter_snps <- function(X, pos, chr = NULL, maf_min = 0.05, callrate_min = 0.8,
                            het_max = 1) {
  callrate <- colMeans(!is.na(X))
  keep <- (ldd_maf(X) >= maf_min) & (callrate >= callrate_min)
  if (is.finite(het_max) && het_max < 1) {
    het <- colMeans(X == 1, na.rm = TRUE)   # heterozygote ratio per marker
    keep <- keep & (het <= het_max)
  }
  keep[is.na(keep)] <- FALSE
  list(X = X[, keep, drop = FALSE], pos = pos[keep],
    chr = if (is.null(chr)) NULL else chr[keep])
}

ldd_compute_ld_pairs <- function(X, pos, max_dist = 5e5, max_snps = 6000L,
                                 window_mode = "distance", max_snp_window = 100L,
                                 min_r2_keep = 0, measure = "r2",
                                 pair_cap = 2e5, seed = 1L) {
  stopifnot(ncol(X) == length(pos))
  if (ncol(X) < 2L) return(data.frame(dist = numeric(0), r2 = numeric(0)))
  ord <- order(pos)
  X <- X[, ord, drop = FALSE]
  pos <- pos[ord]

  # Guard against pathological in-memory sizes: thin evenly if too many SNPs.
  m <- ncol(X)
  if (m > max_snps) {
    idx <- round(seq(1, m, length.out = max_snps))
    X <- X[, idx, drop = FALSE]
    pos <- pos[idx]
    m <- ncol(X)
  }

  dprime <- identical(measure, "dprime")
  Cov <- suppressWarnings(stats::cov(X, use = "pairwise.complete.obs"))
  sdv <- sqrt(diag(Cov))
  R2 <- (Cov / outer(sdv, sdv))^2
  if (dprime) {
    Dm <- Cov / 2                       # composite LD coefficient D_AB
    p <- colMeans(X, na.rm = TRUE) / 2  # frequency of the counted allele
  }

  out <- vector("list", m - 1L)
  for (i in seq_len(m - 1L)) {
    j <- (i + 1L):m
    d <- pos[j] - pos[i]
    keepw <- switch(window_mode,
      snp_count = (j - i) <= max_snp_window,
      both = (d <= max_dist) & ((j - i) <= max_snp_window),
      d <= max_dist)                    # "distance" (default)
    r2ij <- R2[i, j]
    keep <- keepw & is.finite(d) & is.finite(r2ij) & (r2ij >= min_r2_keep)
    if (!any(keep)) next
    if (dprime) {
      Dij <- Dm[i, j]; pj <- p[j]; pi <- p[i]
      Dmax <- ifelse(Dij > 0, pmin(pi * (1 - pj), (1 - pi) * pj),
        pmin(pi * pj, (1 - pi) * (1 - pj)))
      val <- pmin(pmax(abs(Dij) / Dmax, 0), 1)
    } else val <- r2ij
    out[[i]] <- data.frame(dist = d[keep], r2 = val[keep])
  }
  res <- do.call(rbind, out)
  if (is.null(res)) return(data.frame(dist = numeric(0), r2 = numeric(0)))
  res <- res[is.finite(res$dist) & is.finite(res$r2), , drop = FALSE]

  # Reproducible random subsample keeps the distance distribution unbiased.
  if (is.finite(pair_cap) && pair_cap > 0 && nrow(res) > pair_cap) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
    res <- res[sort(sample.int(nrow(res), pair_cap)), , drop = FALSE]
  }
  res
}

ldd_bin_ld <- function(ld, bin_bp = 1000) {
  if (!nrow(ld)) return(data.frame(dist = numeric(0), mean_r2 = numeric(0), n_pairs = integer(0)))
  bin_bp <- max(1, bin_bp)
  bin <- floor(ld$dist / bin_bp) * bin_bp + bin_bp / 2
  agg_mean <- tapply(ld$r2, bin, mean, na.rm = TRUE)
  agg_n <- tapply(ld$r2, bin, length)
  data.frame(dist = as.numeric(names(agg_mean)), mean_r2 = as.numeric(agg_mean),
    n_pairs = as.integer(agg_n), row.names = NULL)[order(as.numeric(names(agg_mean))), ]
}

ldd_auto_bin <- function(max_dist) {
  if (!is.finite(max_dist) || max_dist <= 0) return(1000)
  target <- max_dist / 200          # aim for ~200 bins
  mag <- 10^floor(log10(target))
  step <- c(1, 2, 5, 10)[which.min(abs(c(1, 2, 5, 10) - target / mag))]
  max(1, step * mag)
}

ldd_bin_ld_log <- function(ld, n = 200L) {
  ld <- ld[is.finite(ld$dist) & ld$dist > 0, , drop = FALSE]
  if (!nrow(ld)) return(data.frame(dist = numeric(0), mean_r2 = numeric(0), n_pairs = integer(0)))
  n <- max(2L, as.integer(n))
  lb <- log10(range(ld$dist))
  brk <- seq(lb[1], lb[2], length.out = n + 1L)
  idx <- findInterval(log10(ld$dist), brk, rightmost.closed = TRUE)
  mid <- 10^((brk[-length(brk)] + brk[-1]) / 2)
  m <- tapply(ld$r2, idx, mean, na.rm = TRUE)
  cnt <- tapply(ld$r2, idx, length)
  k <- as.integer(names(m))
  data.frame(dist = mid[k], mean_r2 = as.numeric(m), n_pairs = as.integer(cnt),
    row.names = NULL)[order(mid[k]), ]
}

ldd_resolve_bin <- function(dmax, settings = ldd_default_settings()) {
  mode <- settings$bin_mode %||% "auto"
  if (identical(mode, "fixed_width") && (settings$bin_width_kb %||% 0) > 0)
    return(list(type = "width", width = settings$bin_width_kb * 1000))
  if (identical(mode, "n_bins"))
    return(list(type = "width", width = max(1, dmax / max(1L, settings$n_bins %||% 200L))))
  if (identical(mode, "log"))
    return(list(type = "log", n = settings$n_bins %||% 200L))
  list(type = "width", width = ldd_auto_bin(dmax))
}

ldd_hw_expected_r2 <- function(dist, rho, n) {
  C <- rho * dist
  term1 <- (10 + C) / ((2 + C) * (11 + C))
  term2 <- 1 + ((3 + C) * (12 + 12 * C + C^2)) / (n * (2 + C) * (11 + C))
  term1 * term2
}

ldd_exp_expected_r2 <- function(dist, rho, a) {
  a / (1 + rho * dist)
}

ldd_fit_decay <- function(ld, n_ind, model = "hw") {
  ld <- ld[is.finite(ld$dist) & is.finite(ld$r2) & ld$dist > 0, , drop = FALSE]
  n_ind_value <- suppressWarnings(as.numeric(n_ind))
  n_known <- length(n_ind_value) == 1L && is.finite(n_ind_value) && n_ind_value > 1
  n <- if (n_known) 2 * n_ind_value else NA_real_
  # The finite-sample Hill-Weir correction needs the number of sampled
  # chromosomes. For LD tables without sample metadata, use its large-sample
  # limit instead of fabricating a sample size from the number of LD records.
  n_model <- if (n_known) n else 1e9
  fail <- list(model = model, rho = NA_real_, n = n,
    converged = FALSE, predict = function(d) rep(NA_real_, length(d)), params = NULL)
  if (identical(model, "none") || nrow(ld) < 5L) {
    fail$model <- if (identical(model, "none")) "none" else model
    return(fail)
  }

  rho0 <- 1 / stats::median(ld$dist)
  if (!is.finite(rho0) || rho0 <= 0) return(fail)

  if (identical(model, "hw")) {
    fit <- .ldd_nls_bounded(
      r2 ~ ldd_hw_expected_r2(dist, rho, n_model), ld,
      start = list(rho = rho0),
      lower = c(rho = .LDD_RHO_MIN), upper = c(rho = .LDD_RHO_MAX))
    if (is.null(fit)) return(fail)
    rho <- unname(stats::coef(fit)["rho"])
    predict_fn <- function(d) ldd_hw_expected_r2(d, rho, n_model)
    if (!.ldd_decay_is_valid(rho, predict_fn, ld$dist)) return(fail)
    return(list(model = "hw", rho = rho, n = n, converged = TRUE,
      predict = predict_fn,
      params = list(rho = rho, n = n, n_model = n_model,
        sample_size_assumption = if (n_known) "observed" else "large-sample")))
  }

  # Exponential / Sved. `a` is the distance-zero LD asymptote, so it is bounded
  # to the (0, 1] range that r^2 / D' can actually take.
  a0 <- max(ld$r2, na.rm = TRUE)
  if (!is.finite(a0) || a0 <= 0) return(fail)
  fit <- .ldd_nls_bounded(
    r2 ~ ldd_exp_expected_r2(dist, rho, a), ld,
    start = list(rho = rho0, a = min(a0, 1)),
    lower = c(rho = .LDD_RHO_MIN, a = 1e-6), upper = c(rho = .LDD_RHO_MAX, a = 1))
  if (is.null(fit)) return(fail)
  co <- stats::coef(fit)
  rho <- unname(co["rho"]); a <- unname(co["a"])
  predict_fn <- function(d) ldd_exp_expected_r2(d, rho, a)
  if (!.ldd_decay_is_valid(rho, predict_fn, ld$dist)) return(fail)
  list(model = "exp", rho = rho, n = n, converged = TRUE,
    predict = predict_fn,
    params = list(rho = rho, a = a))
}

.LDD_RHO_MIN <- 1e-12

.LDD_RHO_MAX <- 1

.ldd_nls_bounded <- function(formula, data, start, lower, upper) {
  ctrl <- stats::nls.control(maxiter = 500, warnOnly = TRUE)
  fit <- tryCatch(
    suppressWarnings(stats::nls(formula, data = data, start = start,
      algorithm = "port", lower = lower, upper = upper, control = ctrl)),
    error = function(e) NULL)
  if (is.null(fit))
    fit <- tryCatch(
      suppressWarnings(stats::nls(formula, data = data, start = start, control = ctrl)),
      error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  conv <- tryCatch(isTRUE(fit$convInfo$isConv), error = function(e) FALSE)
  # port reports convergence via stopMessage / convInfo$isConv depending on
  # version; treat a finite objective with finite coefficients as the fallback.
  if (!conv) {
    co <- tryCatch(stats::coef(fit), error = function(e) NULL)
    conv <- !is.null(co) && all(is.finite(co))
  }
  if (!conv) return(NULL)
  co <- stats::coef(fit)
  if (!all(is.finite(co))) return(NULL)
  fit
}

.ldd_decay_is_valid <- function(rho, predict_fn, distances) {
  if (!is.finite(rho) || rho <= .LDD_RHO_MIN) return(FALSE)
  rng <- suppressWarnings(range(distances[is.finite(distances) & distances > 0]))
  if (!all(is.finite(rng)) || rng[2] <= rng[1]) return(FALSE)
  ends <- suppressWarnings(predict_fn(c(max(rng[1], 1), rng[2])))
  if (!all(is.finite(ends))) return(FALSE)
  if (any(ends < -1e-8) || any(ends > 1 + 1e-8)) return(FALSE)
  # Require a visible drop: at least 1% of the starting level across the range.
  ends[1] - ends[2] > max(1e-6, 0.01 * ends[1])
}

ldd_distance_at <- function(fit, target_r2, dist_max = 5e6, n_grid = 20000) {
  if (!isTRUE(fit$converged)) return(NA_real_)
  g <- seq(1, dist_max, length.out = n_grid)
  r2 <- fit$predict(g)
  idx <- which(r2 <= target_r2)[1]
  if (is.na(idx) || idx < 2L) return(NA_real_)
  stats::approx(x = r2[c(idx - 1L, idx)], y = g[c(idx - 1L, idx)], xout = target_r2)$y
}

ldd_half_decay <- function(fit, dist_max = 5e6) {
  if (!isTRUE(fit$converged)) return(NA_real_)
  r2_0 <- fit$predict(1)
  ldd_distance_at(fit, target_r2 = r2_0 / 2, dist_max = dist_max)
}

ldd_metrics <- function(fit, threshold_r2 = 0.1, dist_max = 5e6) {
  list(
    rho = fit$rho,
    r2_0 = if (isTRUE(fit$converged)) fit$predict(1) else NA_real_,
    half_decay = ldd_half_decay(fit, dist_max),
    threshold_r2 = threshold_r2,
    threshold_dist = ldd_distance_at(fit, threshold_r2, dist_max),
    converged = isTRUE(fit$converged))
}

ldd_pair_summary <- function(ld = NULL, binned = NULL) {
  empty <- list(n_records = 0L, n_bins = if (is.null(binned)) 0L else nrow(binned),
    max_dist = NA_real_, mean_ld = NA_real_, median_ld = NA_real_, p95_ld = NA_real_,
    pct_ge_0_2 = NA_real_, pct_ge_0_8 = NA_real_, approximate = FALSE)
  if (!is.null(ld) && nrow(ld)) {
    keep <- is.finite(ld$dist) & is.finite(ld$r2) & ld$dist >= 0
    values <- ld$r2[keep]
    distances <- ld$dist[keep]
    if (!length(values)) return(empty)
    return(list(n_records = length(values), n_bins = if (is.null(binned)) 0L else nrow(binned),
      max_dist = max(distances), mean_ld = mean(values), median_ld = stats::median(values),
      p95_ld = as.numeric(stats::quantile(values, 0.95, names = FALSE)),
      pct_ge_0_2 = 100 * mean(values >= 0.2), pct_ge_0_8 = 100 * mean(values >= 0.8),
      approximate = FALSE))
  }
  if (is.null(binned) || !nrow(binned)) return(empty)
  keep <- is.finite(binned$dist) & is.finite(binned$mean_r2)
  values <- binned$mean_r2[keep]
  distances <- binned$dist[keep]
  weights <- if ("n_pairs" %in% names(binned)) binned$n_pairs[keep] else rep(1, length(values))
  weights[!is.finite(weights) | weights <= 0] <- 1
  if (!length(values)) return(empty)
  weighted_quantile <- function(probability) {
    order_index <- order(values)
    sorted_values <- values[order_index]
    cumulative <- cumsum(weights[order_index]) / sum(weights)
    sorted_values[which(cumulative >= probability)[1]]
  }
  list(n_records = as.integer(round(sum(weights))), n_bins = length(values),
    max_dist = max(distances), mean_ld = stats::weighted.mean(values, weights),
    median_ld = weighted_quantile(0.5), p95_ld = weighted_quantile(0.95),
    pct_ge_0_2 = 100 * sum(weights[values >= 0.2]) / sum(weights),
    pct_ge_0_8 = 100 * sum(weights[values >= 0.8]) / sum(weights), approximate = TRUE)
}

ldd_fit_diagnostics <- function(fit, fit_data, source = "LD records") {
  empty <- list(converged = FALSE, n_fit = 0L, source = source, rmse = NA_real_,
    mae = NA_real_, bias = NA_real_, pseudo_r2 = NA_real_, spearman = NA_real_)
  if (is.null(fit_data) || !nrow(fit_data) || !isTRUE(fit$converged)) return(empty)
  keep <- is.finite(fit_data$dist) & is.finite(fit_data$r2) & fit_data$dist > 0
  observed <- fit_data$r2[keep]
  distance <- fit_data$dist[keep]
  predicted <- tryCatch(fit$predict(distance), error = function(e) rep(NA_real_, length(distance)))
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  distance <- distance[valid]
  if (!length(observed)) return(empty)
  residual <- observed - predicted
  rss <- sum(residual^2)
  tss <- sum((observed - mean(observed))^2)
  list(converged = TRUE, n_fit = length(observed), source = source,
    rmse = sqrt(mean(residual^2)), mae = mean(abs(residual)), bias = mean(residual),
    pseudo_r2 = if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_,
    spearman = suppressWarnings(stats::cor(distance, observed, method = "spearman")))
}

ldd_make_entry <- function(ld = NULL, binned = NULL, n_ind = NULL, settings = ldd_default_settings(),
                           keep_raw = TRUE) {
  dmax <- 5e5
  if (!is.null(ld) && nrow(ld)) dmax <- max(ld$dist)
  else if (!is.null(binned) && nrow(binned)) dmax <- max(binned$dist)
  if (is.null(binned) && !is.null(ld) && nrow(ld)) {
    br <- ldd_resolve_bin(dmax, settings)
    binned <- if (identical(br$type, "log")) ldd_bin_ld_log(ld, br$n) else ldd_bin_ld(ld, br$width)
  }
  use_raw_fit <- !is.null(ld) && nrow(ld) >= 50L
  fit_data <- if (use_raw_fit) ld else if (!is.null(binned) && nrow(binned))
    data.frame(dist = binned$dist, r2 = binned$mean_r2) else
      data.frame(dist = numeric(0), r2 = numeric(0))
  fit_source <- if (use_raw_fit) "pairwise LD records" else "distance-bin means"
  fit <- ldd_fit_decay(fit_data, n_ind, settings$model %||% "hw")
  metrics <- ldd_metrics(fit, settings$threshold_r2 %||% 0.1, dist_max = max(5e6, dmax * 4))
  list(raw = if (keep_raw) ld else NULL, binned = binned, fit = fit,
    metrics = metrics, diagnostics = ldd_fit_diagnostics(fit, fit_data, fit_source),
    pair_summary = ldd_pair_summary(ld, binned), dmax = dmax, n_ind = n_ind)
}

