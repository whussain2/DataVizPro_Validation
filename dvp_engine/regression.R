# ------------------------------------------------------------------------------
# DataVizPro computation engine (extract) - regression
# Source: DataVizPro v1.1.1 (commit af537d9), file modules/regressions/helpers/reg_helpers.R
# Functions are copied VERBATIM from the application; only the subset needed
# to reproduce the manuscript analyses is included. No UI, plotting-server or
# account code is part of this repository.
# Copyright (c) 2026 PlantNura Consulting & Solutions LLP. MIT License (see LICENSE).
# ------------------------------------------------------------------------------

fit_regression_model <- function(data, response, predictors, reg_type, reg_options = list()) {

  tryCatch({
    # Drop rows with NA in the selected columns
    cols <- c(response, predictors)
    data <- data[complete.cases(data[, cols, drop = FALSE]), , drop = FALSE]

    if (nrow(data) < 5) stop("Not enough complete observations (need at least 5).")

    y <- data[[response]]
    x_df <- data[, predictors, drop = FALSE]
    
    # REG-013: refuse clearly wrong response types with a plain message
    is_count_type <- reg_type %in% c("poisson", "quasipoisson", "negbin")
    if (is_count_type) {
      if (!is.numeric(y)) stop(paste0("'", response, "' is not numeric. Count models need whole-number counts (0, 1, 2 ...)."))
      if (any(y < 0, na.rm = TRUE) || any(abs(y - round(y)) > 1e-8, na.rm = TRUE))
        stop(paste0("'", response, "' has decimals or negative values. Count models (Poisson / Negative Binomial) need whole-number counts; use Linear or Gamma for a continuous response."))
    }
    if (reg_type %in% c("multinomial", "ordinal")) {
      n_cls <- length(unique(y))
      if (n_cls > 20) stop(paste0("'", response, "' has ", n_cls, " different values. Multinomial / ordinal models need a response with a few classes (2-20); a continuous response needs Linear or another numeric model."))
      if (n_cls < 2) stop(paste0("'", response, "' has only one class."))
    }
    if (reg_type %in% c("logistic", "probit")) {
      # REG-001: a text (Yes/No) response is turned into a 2-level factor;
      # the second level is the modelled "success"
      if (!is.numeric(y)) {
        lv <- unique(as.character(y[!is.na(y)]))
        if (length(lv) != 2) stop(paste0("'", response, "' has ", length(lv), " classes. Logistic / probit models need exactly 2 (e.g. Yes/No or 0/1)."))
        y_fac <- factor(as.character(y), levels = sort(lv))
        data[[response]] <- y_fac; y <- y_fac
      } else if (!all(y %in% c(0, 1))) {
        stop(paste0("'", response, "' must be 0/1 (or two text classes) for a logistic / probit model."))
      }
    }
    # a model with as many parameters as rows cannot be estimated
    n_par <- sum(vapply(predictors, function(v) if (is.numeric(data[[v]])) 1L else max(1L, length(unique(data[[v]])) - 1L), integer(1))) + 1L
    if (reg_type %in% c("linear", "polynomial", "robust", "logistic", "probit", "poisson", "quasipoisson", "negbin", "gamma", "beta") && n_par >= nrow(data)) {
      stop(paste0("The model has ", n_par, " parameters for ", nrow(data), " rows (a predictor with very many categories?). Choose fewer predictors or a predictor with fewer groups."))
    }
    if (reg_type == "polynomial") {
      deg <- suppressWarnings(as.integer(reg_options$poly_degree))
      if (length(deg) != 1 || is.na(deg) || deg < 2 || deg > 10) stop("Polynomial degree must be a whole number from 2 to 10.")  # REG-021
    }

    switch(reg_type,

      # ---- Linear ----
      "linear" = {
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- lm(as.formula(formula_str), data = data)
        list(model = model, type = "linear", data = data, response = response, predictors = predictors)
      },

      # ---- Polynomial ----
      "polynomial" = {
        deg <- as.integer(reg_options$poly_degree %||% 2)
        # Build polynomial terms for each predictor (supports 1 or more)
        poly_terms <- vapply(predictors, function(pred) {
          paste0("poly(", pred, ", ", deg, ", raw = TRUE)")
        }, character(1))
        formula_str <- paste(response, "~", paste(poly_terms, collapse = " + "))
        model <- lm(as.formula(formula_str), data = data)
        list(model = model, type = "polynomial", degree = deg, data = data, response = response, predictors = predictors)
      },

      # ---- Logistic (Binary) ----
      "logistic" = {
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- glm(as.formula(formula_str), data = data, family = binomial(link = "logit"))
        list(model = model, type = "logistic", data = data, response = response, predictors = predictors)
      },

      # ---- Poisson ----
      "poisson" = {
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- glm(as.formula(formula_str), data = data, family = poisson(link = "log"))
        list(model = model, type = "poisson", data = data, response = response, predictors = predictors)
      },

      # ---- Negative Binomial ----
      "negbin" = {
        if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required for Negative Binomial regression.")
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        nb_warn <- character(0)
        model <- withCallingHandlers(MASS::glm.nb(as.formula(formula_str), data = data),
                                     warning = function(w) { nb_warn <<- c(nb_warn, conditionMessage(w)); invokeRestart("muffleWarning") })
        if (length(nb_warn) > 0) {  # REG-014
          dom <- shiny::getDefaultReactiveDomain()
          if (!is.null(dom)) shiny::showNotification(
            paste0("\u26a0 Negative Binomial: ", paste(unique(nb_warn), collapse = "; "),
                   ". The data may not be overdispersed - compare with the Poisson model."),
            type = "warning", duration = 10, session = dom)
        }
        list(model = model, type = "negbin", data = data, response = response, predictors = predictors)
      },

      # ---- Ridge ----
      "ridge" = {
        if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required for Ridge regression.")
        x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
        y_vec <- as.numeric(y)
        model <- glmnet::glmnet(x_mat, y_vec, alpha = 0, standardize = TRUE)
        cv_model <- glmnet::cv.glmnet(x_mat, y_vec, alpha = 0, nfolds = 5)
        list(model = model, cv_model = cv_model, type = "ridge",
             x_mat = x_mat, y_vec = y_vec, data = data, response = response, predictors = predictors)
      },

      # ---- Lasso ----
      "lasso" = {
        if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required for Lasso regression.")
        x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
        y_vec <- as.numeric(y)
        model <- glmnet::glmnet(x_mat, y_vec, alpha = 1, standardize = TRUE)
        cv_model <- glmnet::cv.glmnet(x_mat, y_vec, alpha = 1, nfolds = 5)
        list(model = model, cv_model = cv_model, type = "lasso",
             x_mat = x_mat, y_vec = y_vec, data = data, response = response, predictors = predictors)
      },

      # ---- Elastic Net ----
      "elasticnet" = {
        if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required for Elastic Net regression.")
        alpha_val <- as.numeric(reg_options$alpha %||% 0.5)
        x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
        y_vec <- as.numeric(y)
        model <- glmnet::glmnet(x_mat, y_vec, alpha = alpha_val, standardize = TRUE)
        cv_model <- glmnet::cv.glmnet(x_mat, y_vec, alpha = alpha_val, nfolds = 5)
        list(model = model, cv_model = cv_model, type = "elasticnet", alpha = alpha_val,
             x_mat = x_mat, y_vec = y_vec, data = data, response = response, predictors = predictors)
      },

      # ---- Quantile ----
      "quantile" = {
        if (!requireNamespace("quantreg", quietly = TRUE)) stop("Package 'quantreg' required for Quantile regression.")
        tau_val <- as.numeric(reg_options$tau %||% 0.5)
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- quantreg::rq(as.formula(formula_str), tau = tau_val, data = data)
        list(model = model, type = "quantile", tau = tau_val, data = data, response = response, predictors = predictors)
      },

      # ---- Robust ----
      "robust" = {
        if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required for Robust regression.")
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- MASS::rlm(as.formula(formula_str), data = data)
        list(model = model, type = "robust", data = data, response = response, predictors = predictors)
      },

      # ---- GAM ----
      "gam" = {
        if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' required for GAM.")
        smooth_terms <- paste(paste0("s(", predictors, ")"), collapse = " + ")
        formula_str <- paste(response, "~", smooth_terms)
        model <- mgcv::gam(as.formula(formula_str), data = data)
        list(model = model, type = "gam", data = data, response = response, predictors = predictors)
      },

      # ---- Gamma (positive continuous, right-skewed) ----
      "gamma" = {
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- glm(as.formula(formula_str), data = data, family = Gamma(link = "log"))
        list(model = model, type = "gamma", data = data, response = response, predictors = predictors)
      },

      # ---- Beta (response in 0-1, proportions / rates) ----
      "beta" = {
        if (!requireNamespace("betareg", quietly = TRUE)) stop("Package 'betareg' required for Beta regression.")
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- betareg::betareg(as.formula(formula_str), data = data)
        list(model = model, type = "beta", data = data, response = response, predictors = predictors)
      },

      # ---- Quasi-Poisson (overdispersed count data) ----
      "quasipoisson" = {
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- glm(as.formula(formula_str), data = data, family = quasipoisson(link = "log"))
        list(model = model, type = "quasipoisson", data = data, response = response, predictors = predictors)
      },

      # ---- Probit (binary response, probit link) ----
      "probit" = {
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        model <- glm(as.formula(formula_str), data = data, family = binomial(link = "probit"))
        list(model = model, type = "probit", data = data, response = response, predictors = predictors)
      },

      # ---- Multinomial Logistic (nominal multi-class response) ----
      "multinomial" = {
        if (!requireNamespace("nnet", quietly = TRUE)) stop("Package 'nnet' required for Multinomial Logistic regression.")
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        data[[response]] <- as.factor(data[[response]])
        model <- nnet::multinom(as.formula(formula_str), data = data, trace = FALSE)
        list(model = model, type = "multinomial", data = data, response = response, predictors = predictors)
      },

      # ---- Ordinal Logistic (ordered multi-class response) ----
      "ordinal" = {
        if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required for Ordinal Logistic regression.")
        formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
        # REG-002: level order. A factor keeps its levels; known ordinal words
        # (Low/Medium/High, Early/Medium/Late, ...) and numbers are put in
        # their natural order; otherwise alphabetical as before.
        data[[response]] <- reg_ordered_response(data[[response]])
        model <- MASS::polr(as.formula(formula_str), data = data, Hess = TRUE)
        list(model = model, type = "ordinal", data = data, response = response, predictors = predictors,
             response_levels = levels(data[[response]]))
      },

      stop(paste("Unknown regression type:", reg_type))
    )

  }, error = function(e) {
    list(error = e$message, type = reg_type)
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

reg_ordered_response <- function(y) {
  if (is.factor(y)) return(factor(as.character(y), levels = levels(droplevels(y)), ordered = TRUE))
  lv <- unique(as.character(y[!is.na(y)]))
  known <- list(
    c("very low", "low", "medium", "moderate", "high", "very high"),
    c("none", "low", "medium", "high"),
    c("early", "medium", "mid", "intermediate", "late"),
    c("small", "medium", "large"),
    c("short", "medium", "intermediate", "tall"),
    c("poor", "fair", "good", "very good", "excellent"),
    c("resistant", "moderately resistant", "moderately susceptible", "susceptible", "highly susceptible"),
    c("never", "rarely", "sometimes", "often", "always"),
    c("strongly disagree", "disagree", "neutral", "agree", "strongly agree")
  )
  low <- tolower(trimws(lv))
  for (k in known) {
    if (all(low %in% k)) return(factor(as.character(y), levels = lv[order(match(low, k))], ordered = TRUE))
  }
  nums <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", lv)))
  if (!anyNA(nums) && any(grepl("[0-9]", lv))) return(factor(as.character(y), levels = lv[order(nums)], ordered = TRUE))
  factor(as.character(y), levels = sort(lv), ordered = TRUE)
}

