# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - fieldspatial
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/heatmaps/helpers/fieldspatial_helpers.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

fs_nz <- function(x, default) {
  if (is.null(x) || length(x) == 0) return(default)
  if (length(x) == 1 && (is.na(x) || (is.character(x) && !nzchar(x)))) return(default)
  x
}

fs_message_plot <- function(msg, subtitle = NULL) {
  df <- data.frame(x = 1, y = 1)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::annotate("text", x = 1, y = 1.06, label = msg,
                      size = 5, fontface = "bold", colour = "#34495e") +
    ggplot2::theme_void() +
    ggplot2::xlim(0, 2) + ggplot2::ylim(0, 2)
  if (!is.null(subtitle) && nzchar(subtitle)) {
    p <- p + ggplot2::annotate("text", x = 1, y = 0.9, label = subtitle,
                               size = 3.6, colour = "#7f8c8d")
  }
  p
}

fs_default_settings <- function() {
  list(
    view_type = "raw",            # raw | smooth | surface3d | trend | autocorr | temporal

    # Appearance
    resolution      = 120,        # smooth / 3D interpolation grid size
    show_contours   = FALSE,
    contour_bins    = 10,
    show_points     = FALSE,
    point_size      = 1.2,
    show_borders    = TRUE,
    border_color    = "white",
    border_size     = 0.3,
    equal_aspect    = TRUE,
    cellwidth       = NA_real_,   # px per cell (NA/0 = auto-fill plot area)
    cellheight      = NA_real_,   # px per cell (NA/0 = auto-fill plot area)
    fill_alpha      = 1,          # tile fill opacity (0-1)
    angle_theta     = 35,         # 3D azimuth
    angle_phi       = 25,         # 3D elevation
    show_3d_contours = TRUE,
    surface3d_expand = 0.6,       # 3D vertical exaggeration (persp expand)
    surface3d_shade  = 0.3,       # 3D surface shading (persp shade; 0 = flat colour)
    surface3d_box    = TRUE,      # 3D bounding box + axis ticks
    surface3d_tick_size  = 0.8,   # 3D tick-label cex
    surface3d_tick_color = "#2c3e50",
    surface3d_tick_face  = "plain",
    surface3d_nticks     = 5,     # 3D axis tick count
    # (3D multi-panel strips / spacing use stripText* / facetPanelSpacing below - AWS-006)
    # Per-axis 3D tick label keys (X/Y share xTickSize/yTickSize; Z needs its own)
    surface3d_z_tick_size   = 10,      # Z tick label size (pt; cex = size/10)
    surface3d_z_tick_color  = "#2c3e50",
    surface3d_z_tick_face   = "plain",
    surface3d_hide_z_ticks  = FALSE,
    surface3d_show_x_tick_marks = TRUE,
    surface3d_show_y_tick_marks = TRUE,
    surface3d_show_z_tick_marks = TRUE,
    surface3d_tick_direction = "out", # out | in | cross
    surface3d_tick_length    = 0.025,   # fraction of projected box extent
    surface3d_tick_mark_color = "#2c3e50",
    surface3d_tick_mark_lty  = "solid",
    surface3d_tick_mark_lwd  = 0.8,
    surface3d_axis_lwd       = 0.8,     # bounding-box line width
    # Axis label keys (also used by ggplot path; defaults used when not yet applied)
    x_label_size  = 10, x_label_color = "#2c3e50",
    x_label_face  = "plain", x_label_font = "sans",
    y_label_size  = 10, y_label_color = "#2c3e50",
    y_label_face  = "plain", y_label_font = "sans",
    temporal_frame_ms = 800,
    temporal_transition_ms = 300,

    # Colours
    palette       = "viridis",
    reverse_colors = FALSE,
    na_color      = "grey92",
    manual_limits = FALSE,
    limit_min     = NA_real_,
    limit_max     = NA_real_,

    # Labels (math-aware via parse_label_safe)
    titleLabel    = "",
    subtitleLabel = "",
    subtitle_size  = 12,
    subtitle_color = "#666666",
    subtitle_font  = "sans",
    subtitle_face  = "plain",
    xLabel        = "",
    yLabel        = "",
    zLabel        = "",
    legendTitle   = "",
    hide_x_label  = FALSE,
    hide_y_label  = FALSE,
    hide_z_label  = FALSE,
    z_label_size  = 10,
    z_label_color = "#2c3e50",
    z_label_face  = "plain",
    z_label_font  = "sans",
    # Axis label position offsets (fraction of bounding-box; 0 = at axis, larger = further)
    x_label_offset = 0.14,
    y_label_offset = 0.14,
    z_label_offset = 0.10,

    # Axis tick labels (Axis Text button — tick mark values)
    xTickSize        = 10,
    xTickColor       = "#2c3e50",
    xTickFace        = "plain",
    xAngle           = 0,
    xAxisMargin      = 5.5,
    hideXTickLabels  = FALSE,
    yTickSize        = 10,
    yTickColor       = "#2c3e50",
    yTickFace        = "plain",
    yAngle           = 0,
    yAxisMargin      = 5.5,
    hideYTickLabels  = FALSE,
    # Text annotations (free-floating)
    fs_annotations = list(),
    # Show values in cells
    show_values    = FALSE,
    value_decimals = 2,
    value_size     = 3,
    value_color    = "#000000",
    xTickColorLine   = "#000000",
    xTickLinetype    = "solid",
    xTickLinewidth   = 0.8,
    yTickLength      = 0.5,
    yTickColorLine   = "#000000",
    yTickLinetype    = "solid",
    yTickLinewidth   = 0.8,
    show_axis_ticks  = FALSE,

    # Theme
    theme_style = "minimal",
    base_size   = 13,
    panel_bg    = "#ffffff",
    plot_bg     = "#ffffff",
    show_grid   = FALSE,
    grid_major  = FALSE,
    grid_minor  = FALSE,
    grid_remove_x = FALSE,
    grid_remove_y = FALSE,

    # Title styling (matches General Heatmap)
    title_size  = 14,
    title_color = "#000000",
    title_font  = "sans",
    title_face  = "bold",
    title_hjust = 0.5,

    # Legend
    show_legend     = TRUE,
    legend_position = "right",
    legend_key_size = 1,

    # Multi-panel / faceting
    facet_by      = character(0),
    facet_nrow    = 0,
    facet_ncol    = 0,            # 0 = auto
    # Strip labels
    stripTextSize    = 10,
    stripTextFace    = "bold",
    stripTextColor   = "#000000",
    stripBgColor     = "#e8e8e8",
    stripBorderColor = "#cccccc",
    stripBorderSize  = 0.5,
    # Panel spacing
    facetPanelSpacing = 4,        # mm between panels
    facet_scales  = "fixed",
    compare_factor = "",
    compare_level_a = "",
    compare_level_b = "",
    compare_active = FALSE,

    # Highlight design-factor outlines (desplot out1 / out2 style)
    highlight_enable  = FALSE,
    highlight_factor1 = "",       # bold outline factor key (e.g. block)
    highlight_factor2 = "",       # thin outline factor key (e.g. rep)
    highlight_color1  = "#111111",
    highlight_color2  = "#d62728",
    highlight_size1   = 1.1,
    highlight_size2   = 0.6,
    highlight_label1  = TRUE,     # show group-level names for primary factor
    highlight_label2  = FALSE,    # show group-level names for secondary factor
    highlight_label_size = 3.5,   # group-name text size (ggplot mm)

    # Spatial model (trend adjustment)
    trend_display = "fitted",     # fitted | residual | adjusted
    trend_k       = 0,            # 0 = auto basis dimension
    geno_random   = FALSE,        # treatment fixed (compare entries) or random (BLUP + H2)
    blocks_effect = "random",     # random | fixed
    blocks_nested = FALSE,        # block WITHIN the preceding factor (alpha-lattice)
    spatial_model = "spline",     # spline | none
    diag_scope    = "design",     # raw | design | full  -- what diagnostics describe

    # Autocorrelation
    autocorr_metric  = "moran",   # moran | lisa | semivariogram
    moran_k          = 4,
    moran_style      = "W",
    lisa_sig         = 0.05,
    lisa_padjust     = "fdr",     # multiple-testing control over the n local tests
    variogram_model  = "Sph",
    variogram_directional = FALSE,
    correlogram_lags = 6
  )
}

fs_build_field <- function(data, x_col, y_col, value_col,
                           factor_cols = character(0), time_col = NULL,
                           geno_col = NULL, block_cols = NULL, unit_cols = NULL,
                           x_step = 1, y_step = 1) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    return(list(ok = FALSE, msg = "No data available."))
  }
  need <- c(x_col, y_col, value_col)
  if (any(is.null(need)) || any(!nzchar(need)) || !all(need %in% names(data))) {
    return(list(ok = FALSE, msg = "Please map Column (X), Row (Y) and Response."))
  }

  if (identical(x_col, y_col))   # HMS-008
    return(list(ok = FALSE, msg = "Column (X) and Row (Y) must be different columns."))
  for (st in list(list(v = x_step, lab = "Plot width (X step)"), list(v = y_step, lab = "Plot length (Y step)"))) {
    v <- suppressWarnings(as.numeric(st$v))
    if (length(v) == 1 && !is.na(v) && (!is.finite(v) || v <= 0))   # HMS-026
      return(list(ok = FALSE, msg = sprintf("%s must be a positive number.", st$lab)))
  }
  df <- data
  df$.x <- suppressWarnings(as.numeric(as.character(df[[x_col]])))
  df$.y <- suppressWarnings(as.numeric(as.character(df[[y_col]])))
  df$.z <- suppressWarnings(as.numeric(as.character(df[[value_col]])))

  ok_rows <- is.finite(df$.x) & is.finite(df$.y)
  df <- df[ok_rows, , drop = FALSE]
  if (nrow(df) == 0) {
    return(list(ok = FALSE,
                msg = "Column (X) and Row (Y) must be numeric coordinates."))
  }
  if (!any(is.finite(df$.z))) {
    return(list(ok = FALSE, msg = "The response column contains no numeric values."))
  }

  # Carry candidate design factors as `.f_<key>` columns. Keys are made
  # syntactically safe so they can be used in facet formulas; the original
  # column name is kept as the display label (factors[[key]] = column name).
  factors <- list()
  factor_cols <- factor_cols[!is.na(factor_cols) & nzchar(factor_cols)]
  factor_cols <- setdiff(unique(factor_cols), c(x_col, y_col, value_col))
  factor_cols <- factor_cols[factor_cols %in% names(data)]
  if (length(factor_cols)) {
    keys <- make.unique(make.names(factor_cols), sep = "_")
    for (i in seq_along(factor_cols)) {
      col <- factor_cols[i]; key <- keys[i]
      df[[paste0(".f_", key)]] <- as.factor(as.character(df[[col]]))
      factors[[key]] <- col
    }
  }

  clean <- function(v) {
    v <- v[!is.na(v) & nzchar(v)]
    v <- setdiff(unique(v), c(x_col, y_col, value_col))
    v[v %in% names(df)]
  }

  geno_col <- clean(geno_col)[1]
  if (!length(geno_col) || is.na(geno_col)) geno_col <- NULL
  if (!is.null(geno_col)) df$.g <- as.factor(as.character(df[[geno_col]]))

  block_cols <- clean(block_cols)
  blocks <- list()
  if (length(block_cols)) {
    bkeys <- make.unique(make.names(block_cols), sep = "_")
    for (i in seq_along(block_cols)) {
      df[[paste0(".b_", bkeys[i])]] <- as.factor(as.character(df[[block_cols[i]]]))
      blocks[[bkeys[i]]] <- block_cols[i]
    }
  }

  unit_cols <- clean(unit_cols)
  unit_key <- NULL; unit_label <- NULL
  if (length(unit_cols)) {
    lab <- do.call(paste, c(lapply(unit_cols, function(cn) as.character(df[[cn]])), sep = " | "))
    df$.u <- factor(lab, levels = unique(lab))

    # Register the field identity as a first-class factor so the existing facet,
    # outline and label machinery can panel by it without special-casing.
    # A single-column identity is usually already a candidate factor -- reuse
    # that key rather than adding a duplicate entry to every picker. A composite
    # identity (site + year) has no single source column, so it gets its own
    # synthetic key, made unique against the real factors.
    unit_label <- paste(unit_cols, collapse = " + ")
    existing   <- names(factors)[vapply(factors, identical, logical(1), unit_cols[1])]
    if (length(unit_cols) == 1 && length(existing)) {
      unit_key <- existing[1]
    } else {
      unit_key <- make.names(if (length(unit_cols) == 1) unit_cols[1] else "Field")
      while (unit_key %in% names(factors)) unit_key <- paste0(unit_key, ".")
      df[[paste0(".f_", unit_key)]] <- df$.u
      factors[[unit_key]] <- unit_label
    }
  } else {
    df$.u <- factor("All data")
  }

  time_present <- !is.null(time_col) && nzchar(time_col) && time_col %in% names(data)
  if (time_present) df$.t <- as.character(df[[time_col]])

  x_step <- suppressWarnings(as.numeric(x_step)); if (!isTRUE(is.finite(x_step)) || x_step <= 0) x_step <- 1
  y_step <- suppressWarnings(as.numeric(y_step)); if (!isTRUE(is.finite(y_step)) || y_step <= 0) y_step <- 1

  out <- list(
    ok = TRUE, df = df,
    x_col = x_col, y_col = y_col, value_col = value_col,
    geno_col = geno_col, block_cols = blocks, unit_cols = unit_cols,
    unit_key = unit_key, unit_label = unit_label,
    x_step = x_step, y_step = y_step,
    factors = factors,
    time_col = if (time_present) time_col else NULL,
    n = nrow(df),
    n_plots = nrow(unique(df[c(".u", ".x", ".y")])),   # HMS-033: positions, not records
    nx = length(unique(df$.x)),
    ny = length(unique(df$.y))
  )
  out$unit_levels <- levels(df$.u)
  out$dup <- fs_duplicate_report(out)
  out$reps_per_cell <- out$dup$max_per_cell
  out
}

fs_duplicate_report <- function(field) {
  df <- field$df
  key <- paste(as.character(df$.u), df$.x, df$.y, sep = "\r")
  tb  <- table(key)
  worst <- if (length(tb)) max(tb) else 1L
  list(
    max_per_cell = as.integer(worst),
    n_cells      = length(tb),
    n_duplicated = as.integer(sum(tb > 1L)),
    clean        = worst <= 1L
  )
}

fs_analysis_points <- function(field) {
  df <- field$df
  ok <- is.finite(df$.z)
  agg <- stats::aggregate(list(z = df$.z[ok]),
                          by = list(x = df$.x[ok], y = df$.y[ok]),
                          FUN = function(v) mean(v, na.rm = TRUE))
  agg <- agg[is.finite(agg$z), , drop = FALSE]
  attr(agg, "reps_per_cell") <- field$reps_per_cell %||% 1L
  agg
}

fs_metric_coords <- function(pts, field) {
  cbind(x = pts$x * (field$x_step %||% 1), y = pts$y * (field$y_step %||% 1))
}

fs_metric_unit <- function(field) {
  if ((field$x_step %||% 1) == 1 && (field$y_step %||% 1) == 1) "plot units" else "m"
}

fs_model_terms <- function(field, ps) {
  df <- field$df
  has_geno <- !is.null(field$geno_col) && !is.null(df$.g) && nlevels(droplevels(df$.g)) >= 2
  bkeys <- names(field$block_cols %||% list())
  bkeys <- bkeys[vapply(bkeys, function(k) {
    v <- df[[paste0(".b_", k)]]
    !is.null(v) && nlevels(droplevels(v)) >= 2
  }, logical(1))]
  list(has_geno = has_geno, block_keys = bkeys)
}

fs_fit_model <- function(field, ps) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(list(ok = FALSE, msg = "The spatial model needs the 'mgcv' package."))
  }
  df <- field$df
  tm <- fs_model_terms(field, ps)

  dat <- data.frame(x = df$.x, y = df$.y, z = df$.z)
  if (tm$has_geno) dat$geno <- droplevels(df$.g)

  # Nested blocking (e.g. block WITHIN rep) is expressed by composing the level
  # labels cumulatively, so an alpha-lattice gets rep + rep:block as intended.
  nested <- isTRUE(ps$blocks_nested)
  blk_names <- character(0)
  if (length(tm$block_keys)) {
    acc <- NULL
    for (i in seq_along(tm$block_keys)) {
      k <- tm$block_keys[i]
      v <- as.character(df[[paste0(".b_", k)]])
      lab <- if (nested && !is.null(acc)) paste(acc, v, sep = ":") else v
      acc <- lab
      nm <- paste0("blk", i)
      dat[[nm]] <- factor(lab)
      blk_names <- c(blk_names, nm)
    }
    blk_names <- blk_names[vapply(blk_names, function(n) nlevels(dat[[n]]) >= 2, logical(1))]
  }

  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < 6) return(list(ok = FALSE, msg = "Too few plots to fit a model."))
  nx <- length(unique(dat$x)); ny <- length(unique(dat$y))

  geno_random  <- isTRUE(ps$geno_random)
  blocks_random <- !identical(ps$blocks_effect, "fixed")
  want_spatial <- !identical(ps$spatial_model, "none") && nx >= 3 && ny >= 3

  k_user <- fs_nz(ps$trend_k, 0)
  kx <- if (k_user > 0) min(k_user, nx) else max(3, min(10, nx - 1))
  ky <- if (k_user > 0) min(k_user, ny) else max(3, min(10, ny - 1))

  design_terms <- c(
    if (tm$has_geno) { if (geno_random) "s(geno, bs = 're')" else "geno" },
    vapply(blk_names, function(n)
      if (blocks_random) sprintf("s(%s, bs = 're')", n) else n, character(1), USE.NAMES = FALSE)
  )
  sp_term <- function(kx, ky) sprintf("te(x, y, bs = 'ps', k = c(%d, %d))", kx, ky)
  build <- function(terms) stats::as.formula(
    paste("z ~", if (length(terms)) paste(terms, collapse = " + ") else "1"))

  gam_try <- function(terms) tryCatch(
    suppressWarnings(mgcv::gam(build(terms), data = dat, method = "REML")),
    error = function(e) NULL)

  # Step the basis down before dropping design effects: the design is what the
  # user declared, the smooth size is only a tuning choice.
  fit <- NULL; used_design <- design_terms; used_k <- c(kx, ky)
  if (want_spatial) {
    for (kk in list(c(kx, ky), c(max(3, kx - 2), max(3, ky - 2)), c(3, 3))) {
      fit <- gam_try(c(design_terms, sp_term(kk[1], kk[2])))
      if (!is.null(fit)) { used_k <- kk; break }
    }
    if (is.null(fit)) {
      fit <- gam_try(c(design_terms, sp_term(3, 3)))
      if (is.null(fit)) { fit <- gam_try(sp_term(3, 3)); used_design <- character(0) }
    }
  } else {
    fit <- gam_try(design_terms)
  }
  if (is.null(fit)) return(list(ok = FALSE, msg = "The model failed to converge."))

  base_fit <- gam_try(used_design)   # same design, no spatial term

  terms_mat <- tryCatch(stats::predict(fit, type = "terms"), error = function(e) NULL)
  col_of <- function(pat) {
    if (is.null(terms_mat)) return(NULL)
    j <- grep(pat, colnames(terms_mat))
    if (!length(j)) NULL else as.numeric(terms_mat[, j[1]])
  }
  spatial <- if (want_spatial) col_of("^te\\(") else NULL
  if (is.null(spatial)) spatial <- rep(0, nrow(dat))

  block_eff <- rep(0, nrow(dat))
  for (n in blk_names) {
    v <- col_of(sprintf("^s?\\(?%s", n))
    if (!is.null(v)) block_eff <- block_eff + v
  }

  full_fit <- as.numeric(stats::predict(fit))
  list(
    ok = TRUE, fit = fit, base_fit = base_fit, data = dat,
    spatial = spatial, block_eff = block_eff, fitted = full_fit,
    resid = dat$z - full_fit,
    has_geno = tm$has_geno, geno_random = geno_random,
    blocks_random = blocks_random, block_names = blk_names,
    block_labels = unname(unlist(field$block_cols %||% list()))[seq_along(blk_names)],
    spatial_used = want_spatial, k = used_k,
    design_dropped = length(used_design) < length(design_terms),
    n = nrow(dat)
  )
}

fs_diag_points <- function(field, ps) {
  scope <- fs_nz(ps$diag_scope, "design")
  if (identical(scope, "raw")) {
    pts <- fs_analysis_points(field)
    attr(pts, "scope_label") <- "raw values"
    return(pts)
  }
  ps2 <- ps
  if (identical(scope, "design")) ps2$spatial_model <- "none"
  m <- fs_fit_model(field, ps2)
  if (!isTRUE(m$ok) || (!m$has_geno && !length(m$block_names) && identical(scope, "design"))) {
    pts <- fs_analysis_points(field)
    attr(pts, "scope_label") <- if (isTRUE(m$ok))
      "raw values (no design effects were mapped)" else "raw values (model unavailable)"
    return(pts)
  }
  d <- m$data
  agg <- stats::aggregate(list(z = m$resid), by = list(x = d$x, y = d$y), FUN = mean)
  agg <- agg[is.finite(agg$z), , drop = FALSE]
  attr(agg, "reps_per_cell") <- field$reps_per_cell %||% 1L
  attr(agg, "scope_label") <- if (identical(scope, "design"))
    "residuals after design effects" else "residuals after design + spatial"
  agg
}

fs_fit_trend <- function(field, ps) {
  m <- fs_fit_model(field, ps)
  if (!isTRUE(m$ok)) return(list(ok = FALSE, msg = m$msg %||% "Spatial model unavailable."))

  dat <- m$data
  out <- data.frame(
    x        = dat$x,
    y        = dat$y,
    observed = dat$z,
    fitted   = m$spatial,
    residual = m$resid,
    adjusted = dat$z - m$spatial
  )
  out <- stats::aggregate(cbind(observed, fitted, residual, adjusted) ~ x + y,
                          data = out, FUN = mean)

  sm <- summary(m$fit)
  var_z <- stats::var(dat$z, na.rm = TRUE)
  pct_spatial <- if (is.finite(var_z) && var_z > 0)
    100 * stats::var(m$spatial, na.rm = TRUE) / var_z else NA_real_

  parts <- c(
    if (m$has_geno) sprintf("%s (%s)", field$geno_col, if (m$geno_random) "random" else "fixed"),
    if (length(m$block_names)) sprintf("%s (%s)",
      paste(m$block_labels, collapse = if (isTRUE(ps$blocks_nested)) " / " else " + "),
      if (m$blocks_random) "random" else "fixed"),
    if (m$spatial_used) sprintf("te(x, y) P-spline k=%d\u00d7%d", m$k[1], m$k[2])
  )
  engine_lbl <- if (length(parts)) paste(parts, collapse = "  +  ") else "intercept only"

  summ <- list(
    `Deviance explained`     = paste0(round(100 * sm$dev.expl, 1), "%"),
    `Adj. R-squared`         = round(sm$r.sq, 3),
    `Spatial variance share` = if (is.finite(pct_spatial)) paste0(round(pct_spatial, 1), "%") else "\u2014",
    # HMS-030: plots are field positions; records can be several per plot
    `Plots used`             = nrow(unique(dat[intersect(c("x", "y", "u", ".u"), names(dat))])),
    `Records used`           = m$n
  )
  if (!is.null(m$base_fit) && m$spatial_used) {
    d_aic <- AIC(m$base_fit) - AIC(m$fit)
    summ[["Spatial term worth it?"]] <- sprintf("%s (\u0394AIC = %.1f vs no-spatial model)",
      if (is.finite(d_aic) && d_aic > 2) "Yes" else "Marginal", d_aic)
  }
  if (!m$has_geno && !length(m$block_names)) {
    summ[["Note"]] <- "No treatment or blocking mapped \u2014 the smooth can absorb real treatment effects."
  } else if (isTRUE(m$design_dropped)) {
    summ[["Note"]] <- "Some declared effects had to be dropped for the model to fit."
  }

  list(ok = TRUE, data = out, engine = engine_lbl, summary = summ,
       grand = mean(dat$z, na.rm = TRUE), model = m,
       has_design = m$has_geno || length(m$block_names) > 0)
}

fs_morans_i <- function(field, ps) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    return(list(ok = FALSE, msg = "Moran's I needs the 'spdep' package."))
  }
  pts <- fs_diag_points(field, ps)
  if (nrow(pts) < 5) return(list(ok = FALSE, msg = "Too few cells for Moran's I."))
  k <- max(1, min(fs_nz(ps$moran_k, 4), nrow(pts) - 1))
  res <- tryCatch({
    coords <- fs_metric_coords(pts, field)
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
    lw <- spdep::nb2listw(nb, style = fs_nz(ps$moran_style, "W"))
    # Two-sided: negative (dispersed) autocorrelation must not report p = 1.
    mt <- spdep::moran.test(pts$z, lw, zero.policy = TRUE, alternative = "two.sided")
    list(
      I = unname(mt$estimate[1]), E = unname(mt$estimate[2]),
      V = unname(mt$estimate[3]), z = unname(mt$statistic),
      p = mt$p.value, k = k, n = nrow(pts)
    )
  }, error = function(e) list(err = conditionMessage(e)))
  if (!is.null(res$err)) return(list(ok = FALSE, msg = paste("Moran's I failed:", res$err)))

  strength <- if (abs(res$I) < 0.1) "negligible" else if (abs(res$I) < 0.3) "weak" else
    if (abs(res$I) < 0.5) "moderate" else "strong"
  dir <- if (res$I > res$E) "positive (clustered)" else "negative (dispersed)"
  res$interpretation <- sprintf(
    "Moran's I = %.3f (E[I] = %.3f). %s %s autocorrelation; z = %.2f, p = %s.",
    res$I, res$E, strength, dir, res$z, format.pval(res$p, digits = 3, eps = 1e-4))
  res$ok <- TRUE
  res
}

fs_lisa <- function(field, ps) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    return(list(ok = FALSE, msg = "LISA needs the 'spdep' package.",
                plot = fs_message_plot("LISA needs the 'spdep' package")))
  }
  pts <- fs_diag_points(field, ps)
  if (nrow(pts) < 5) {
    return(list(ok = FALSE, msg = "Too few cells for LISA.",
                plot = fs_message_plot("Too few cells for LISA")))
  }
  k <- max(1, min(fs_nz(ps$moran_k, 4), nrow(pts) - 1))
  out <- tryCatch({
    coords <- fs_metric_coords(pts, field)
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
    lw <- spdep::nb2listw(nb, style = fs_nz(ps$moran_style, "W"))
    locm <- spdep::localmoran(pts$z, lw, zero.policy = TRUE)
    zc <- as.numeric(scale(pts$z))
    lagz <- spdep::lag.listw(lw, zc)
    # One local test PER CELL: without an adjustment ~5% of cells are flagged by
    # chance alone, which on a field map reads as real hotspots.
    p_raw  <- locm[, ncol(locm)]
    adj_m  <- fs_nz(ps$lisa_padjust, "fdr")
    pcol   <- if (identical(adj_m, "none")) p_raw else stats::p.adjust(p_raw, method = adj_m)
    sig <- fs_nz(ps$lisa_sig, 0.05)
    cl <- rep("Not significant", nrow(pts))
    cl[zc > 0 & lagz > 0 & pcol <= sig] <- "High-High (hotspot)"
    cl[zc < 0 & lagz < 0 & pcol <= sig] <- "Low-Low (coldspot)"
    cl[zc > 0 & lagz < 0 & pcol <= sig] <- "High-Low (outlier)"
    cl[zc < 0 & lagz > 0 & pcol <= sig] <- "Low-High (outlier)"
    pts$cluster <- factor(cl, levels = c("High-High (hotspot)", "Low-Low (coldspot)",
                                         "High-Low (outlier)", "Low-High (outlier)",
                                         "Not significant"))
    pts$p_raw      <- p_raw
    pts$p_adjusted <- pcol
    pts
  }, error = function(e) NULL)
  if (is.null(out)) {
    return(list(ok = FALSE, msg = "LISA computation failed.",
                plot = fs_message_plot("LISA computation failed")))
  }

  pal <- c("High-High (hotspot)" = "#c0392b", "Low-Low (coldspot)" = "#2471a3",
           "High-Low (outlier)" = "#e67e22", "Low-High (outlier)" = "#27ae60",
           "Not significant" = "grey88")
  xlab <- fs_nz(ps$xLabel, field$x_col)
  ylab_str <- fs_nz(ps$yLabel, field$y_col)
  n_hot  <- sum(out$cluster == "High-High (hotspot)")
  n_cold <- sum(out$cluster == "Low-Low (coldspot)")
  sig_p  <- fs_nz(ps$lisa_sig, 0.05)
  adj_m  <- fs_nz(ps$lisa_padjust, "fdr")
  adj_lbl <- switch(adj_m,
    none = "unadjusted p", bonferroni = "Bonferroni-adjusted p",
    holm = "Holm-adjusted p", "FDR-adjusted p (Benjamini\u2013Hochberg)")
  n_other <- sum(!out$cluster %in% c("High-High (hotspot)", "Low-Low (coldspot)", "Not significant"))
  p <- ggplot2::ggplot(out, ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                         fill = .data[["cluster"]])) +
    ggplot2::geom_tile(colour = "grey92", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = pal, name = "LISA cluster", drop = FALSE) +
    # Grid positions, so label the actual cells rather than fractional ticks.
    ggplot2::scale_x_continuous(expand = c(0.02, 0), breaks = sort(unique(out$x))) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0), breaks = sort(unique(out$y))) +
    ggplot2::labs(
      title    = "LISA Hotspot / Coldspot Map",
      subtitle = sprintf("%s \u2264 %.2g \u2014 %d hotspot, %d coldspot, %d other significant cells",
                         adj_lbl, sig_p, n_hot, n_cold, n_other),
      # An all-grey map is a real result, not a failure — say so on the plot.
      caption  = if (n_hot + n_cold + n_other == 0)
        "No cell reaches significance: this field shows no local clustering strong enough to stand out from its neighbours."
        else NULL,
      x        = xlab,
      y        = ylab_str
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(colour = "grey70", linewidth = 0.3))) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#5d6d7e"),
      plot.caption  = ggplot2::element_text(size = 9, colour = "#7f8c8d", hjust = 0),
      axis.title    = ggplot2::element_text(size = 11, colour = "#2c3e50"),
      axis.text     = ggplot2::element_text(size = 10, colour = "#2c3e50"),
      legend.title  = ggplot2::element_text(size = 10, face = "bold"),
      legend.text   = ggplot2::element_text(size = 9),
      panel.grid    = ggplot2::element_blank(),
      panel.border  = ggplot2::element_rect(colour = "grey70", fill = NA),
      aspect.ratio  = diff(range(out$y)) / diff(range(out$x))
    )
  list(ok = TRUE, plot = p, data = out,
       summary = sprintf("Significant (%s \u2264 %.2g): %d hotspot, %d coldspot cells (k = %d neighbours).",
                         adj_lbl, sig_p, n_hot, n_cold, k))
}

fs_semivariogram <- function(field, ps) {
  if (!requireNamespace("gstat", quietly = TRUE) || !requireNamespace("sp", quietly = TRUE)) {
    return(list(ok = FALSE, msg = "Semivariogram needs the 'gstat' and 'sp' packages.",
                plot = fs_message_plot("Semivariogram needs 'gstat' + 'sp'")))
  }
  pts <- fs_diag_points(field, ps)
  if (nrow(pts) < 8) {
    return(list(ok = FALSE, msg = "Too few cells for a semivariogram.",
                plot = fs_message_plot("Too few cells for a semivariogram")))
  }
  res <- tryCatch({
    # Distances must be on the ground scale, otherwise the range is reported in
    # meaningless index units and a long narrow plot is treated as square.
    sp_df <- data.frame(x = pts$x * (field$x_step %||% 1),
                        y = pts$y * (field$y_step %||% 1),
                        z = pts$z)
    sp::coordinates(sp_df) <- ~ x + y
    directional <- isTRUE(ps$variogram_directional)
    v <- if (directional) {
      gstat::variogram(z ~ 1, data = sp_df, alpha = c(0, 45, 90, 135))
    } else {
      gstat::variogram(z ~ 1, data = sp_df)
    }
    model_type <- fs_nz(ps$variogram_model, "Sph")
    init <- gstat::vgm(psill = stats::var(pts$z) * 0.7, model = model_type,
                       range = max(v$dist) / 2, nugget = stats::var(pts$z) * 0.1)
    fit <- tryCatch(gstat::fit.variogram(v, init), error = function(e) NULL)
    list(v = v, fit = fit, directional = directional)
  }, error = function(e) list(err = conditionMessage(e)))
  if (!is.null(res$err)) {
    return(list(ok = FALSE, msg = paste("Semivariogram failed:", res$err),
                plot = fs_message_plot("Semivariogram failed")))
  }

  v <- res$v
  max_h <- max(v$dist)
  pts_size <- if (nrow(v) <= 15) 3.2 else 2.4
  p <- ggplot2::ggplot(v, ggplot2::aes(x = .data[["dist"]], y = .data[["gamma"]])) +
    ggplot2::geom_point(size = pts_size, colour = "#2c3e50", shape = 16)
  nugget <- sill <- rng <- NA_real_
  eff_rng <- NA_real_
  singular <- FALSE
  extrapolated <- FALSE
  if (!is.null(res$fit) && !res$directional) {
    line <- tryCatch(gstat::variogramLine(res$fit, maxdist = max_h),
                     error = function(e) NULL)
    if (!is.null(line)) {
      nugget <- res$fit$psill[1]
      sill   <- sum(res$fit$psill)
      rng    <- res$fit$range[nrow(res$fit)]
      singular <- isTRUE(attr(res$fit, "singular"))
      mdl    <- as.character(res$fit$model[nrow(res$fit)])
      # Exponential and Gaussian report a range PARAMETER; the distance at which
      # ~95% of the sill is reached is 3a and sqrt(3)a respectively. Showing the
      # bare parameter as "the range" understates the true extent of correlation.
      eff_rng <- switch(mdl, Exp = 3 * rng, Gau = sqrt(3) * rng, rng)
      # The variogram is only computed out to max_h. A range or sill beyond that
      # is an extrapolation: drawing it would stretch the axes until the actual
      # points collapse into a corner and imply a plateau the data never showed.
      extrapolated <- !is.finite(eff_rng) || eff_rng > max_h
      rng_lbl <- if (isTRUE(all.equal(eff_rng, rng))) sprintf("Range = %.3g", rng)
                 else sprintf("Eff. range = %.3g", eff_rng)
      p <- p +
        ggplot2::geom_line(data = line,
                           ggplot2::aes(x = .data[["dist"]], y = .data[["gamma"]]),
                           colour = "#c0392b", linewidth = 1.1) +
        ggplot2::geom_hline(yintercept = nugget, linetype = "dashed", colour = "#7f8c8d", linewidth = 0.6)
      if (!extrapolated) {
        p <- p +
          ggplot2::geom_hline(yintercept = sill, linetype = "dotted", colour = "#7f8c8d", linewidth = 0.6) +
          ggplot2::geom_vline(xintercept = eff_rng, linetype = "dashed", colour = "#2980b9", linewidth = 0.6) +
          ggplot2::annotate("text", x = max_h, y = sill,
            label = sprintf("Sill = %.3g", sill), hjust = 1.05, vjust = -0.5,
            size = 3.3, colour = "#7f8c8d", fontface = "italic") +
          ggplot2::annotate("text", x = eff_rng, y = min(v$gamma),
            label = rng_lbl, hjust = -0.08, vjust = 0.4,
            angle = 90, size = 3.2, colour = "#2980b9", fontface = "italic")
      }
      if (is.finite(nugget) && nugget > 0)
        p <- p + ggplot2::annotate("text", x = 0, y = nugget,
          label = sprintf("Nugget = %.3g", nugget), hjust = -0.08, vjust = -0.4,
          size = 3.2, colour = "#7f8c8d", fontface = "italic")
    }
  }
  if (res$directional) {
    p <- ggplot2::ggplot(v, ggplot2::aes(x = .data[["dist"]], y = .data[["gamma"]],
                                         colour = factor(.data[["dir.hor"]]),
                                         group  = factor(.data[["dir.hor"]]))) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_colour_brewer(palette = "Set1", name = "Direction (°)")
  }
  model_lbl <- if (!is.null(res$fit) && !res$directional) {
    sprintf("Model: %s", fs_nz(ps$variogram_model, "Sph"))
  } else NULL
  p <- p +
    ggplot2::labs(
      title    = "Empirical Semivariogram",
      subtitle = paste(c(if (!is.null(model_lbl)) model_lbl else "Omnidirectional",
                         attr(pts, "scope_label")), collapse = "  \u00b7  "),
      x        = sprintf("Separation distance (%s)", fs_metric_unit(field)),
      y        = "Semivariance  \u03B3(h)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#5d6d7e"),
      axis.title    = ggplot2::element_text(size = 11, colour = "#2c3e50"),
      axis.text     = ggplot2::element_text(size = 10, colour = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border  = ggplot2::element_rect(colour = "grey70", fill = NA),
      legend.title  = ggplot2::element_text(size = 10, face = "bold")
    )

  summary_txt <- if (all(is.finite(c(nugget, sill, rng)))) {
    paste0(
      sprintf("Model %s \u2014 nugget = %.3g, sill = %.3g, effective range = %.3g %s.",
              fs_nz(ps$variogram_model, "Sph"), nugget, sill, eff_rng, fs_metric_unit(field)),
      if (extrapolated)
        sprintf(paste0(" The variogram was only computed out to %.3g %s and has not reached a plateau, ",
                       "so the range and sill are extrapolated beyond the data and are not drawn on the plot ",
                       "\u2014 read them as a lower bound on how far the correlation extends."),
                max_h, fs_metric_unit(field)) else "",
      if (singular) " Warning: the fit is singular, so these parameters are unreliable \u2014 try another model or more plots." else "")
  } else if (res$directional) {
    "Directional semivariogram (0\u00b0/45\u00b0/90\u00b0/135\u00b0)."
  } else {
    "Empirical semivariogram (model fit unavailable)."
  }
  list(ok = TRUE, plot = p, nugget = nugget, sill = sill, range = rng,
       eff_range = eff_rng, singular = singular, summary = summary_txt)
}

