# Internal helpers for storing and displaying the best final model selected by
# tune_MEMWAS(). These helpers intentionally use only base R functionality.

.memwas_deparse_one_line <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = " "))
  if (is.character(x)) return(paste(x, collapse = " "))
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

.memwas_format_number <- function(x, digits = 6L) {
  if (length(x) == 0L) return("NA")
  x <- suppressWarnings(as.numeric(x)[1L])
  if (is.na(x)) return("NA")
  if (!is.finite(x)) return(as.character(x))
  format(signif(x, digits = digits), trim = TRUE)
}

.memwas_format_named_numeric <- function(x, digits = 6L) {
  if (is.null(x) || length(x) == 0L) return("none")
  nms <- names(x)
  x <- as.numeric(x)
  vals <- vapply(x, .memwas_format_number, character(1L), digits = digits)
  if (!is.null(nms) && length(nms) == length(vals) && any(nzchar(nms))) {
    vals <- ifelse(nzchar(nms), paste0(nms, "=", vals), vals)
  }
  paste(vals, collapse = ", ")
}

.memwas_empty_table <- function(...) {
  cols <- list(...)
  out <- as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
  if (length(cols) == 0L) out <- data.frame(stringsAsFactors = FALSE)
  out[0L, , drop = FALSE]
}

.memwas_first_rows <- function(x, n = 10L) {
  if (is.null(x)) return(x)
  n <- as.integer(n[1L])
  if (!is.finite(n) || n < 1L) n <- 10L
  if (is.data.frame(x)) {
    if (nrow(x) == 0L) return(x)
    return(x[seq_len(min(nrow(x), n)), , drop = FALSE])
  }
  if (length(x) == 0L) return(x)
  x[seq_len(min(length(x), n))]
}

.memwas_cat_wrapped <- function(label, text, indent = 2L, exdent = 4L) {
  text <- if (length(text) == 0L || is.na(text[1L]) || !nzchar(text[1L])) "NA" else as.character(text[1L])
  prefix <- paste(rep(" ", indent), collapse = "")
  ex_prefix <- paste(rep(" ", exdent), collapse = "")
  first_prefix <- paste0(prefix, label)
  available <- max(20L, getOption("width", 80L) - nchar(first_prefix))
  wrapped <- strwrap(text, width = available)
  if (length(wrapped) == 0L) wrapped <- ""
  cat(first_prefix, wrapped[1L], "\n", sep = "")
  if (length(wrapped) > 1L) {
    for (i in 2L:length(wrapped)) cat(ex_prefix, wrapped[i], "\n", sep = "")
  }
}

.memwas_fixed_effect_coefficient_table <- function(fit) {
  if (is.null(fit) || is.null(fit$coefficients)) {
    return(.memwas_empty_table(term = character(0L), model_term = character(0L),
                               dot_target = character(0L), estimate = numeric(0L),
                               std_error = numeric(0L), dot_threshold = numeric(0L),
                               minimal_important_value = numeric(0L),
                               threshold = numeric(0L), null_value = numeric(0L),
                               test_threshold = numeric(0L), null_boundary = numeric(0L),
                               alternative = character(0L),
                               direction = character(0L), alpha = numeric(0L),
                               hypothesis = character(0L), test = character(0L),
                               statistic = numeric(0L), p_value = numeric(0L),
                               significant = logical(0L), DOT_label = character(0L),
                               dot_matched_by = character(0L),
                               penalty_factor = numeric(0L), penalized = logical(0L),
                               selected = logical(0L)))
  }
  beta <- fit$coefficients
  if (!is.null(fit$coefficient_table) && is.data.frame(fit$coefficient_table)) {
    tab <- fit$coefficient_table
  } else {
    tab <- data.frame(term = names(beta) %||% paste0("beta", seq_along(beta)),
                      estimate = as.numeric(beta),
                      stringsAsFactors = FALSE,
                      check.names = FALSE)
  }
  if (!"term" %in% names(tab)) tab$term <- names(beta) %||% paste0("beta", seq_along(beta))
  if (!"estimate" %in% names(tab)) tab$estimate <- as.numeric(beta)
  tab$penalty_factor <- ifelse(tab$term == "(Intercept)", 0, 1)
  tab$penalized <- tab$penalty_factor > 0
  tab$selected <- abs(as.numeric(tab$estimate)) > 1e-8
  row.names(tab) <- NULL
  tab
}

.memwas_fixed_formula_with_coefficients <- function(fit, digits = 6L) {
  if (is.null(fit) || is.null(fit$coefficients)) return(NA_character_)
  beta <- as.numeric(fit$coefficients)
  terms <- names(fit$coefficients)
  if (is.null(terms)) terms <- paste0("beta", seq_along(beta))
  response <- "eta"
  form <- fit$formal_formula %||% fit$formula
  if (inherits(form, "formula") && length(form) >= 2L) {
    response <- paste(deparse(form[[2L]], width.cutoff = 500L), collapse = " ")
  }
  family <- .memwas_normalize_family(fit$family %||% "gaussian")
  lhs <- switch(family,
                gaussian = paste0("E(", response, ")"),
                binomial = paste0("logit(P(", response, " = 1))"),
                poisson = paste0("log(E(", response, "))"),
                negative_binomial = paste0("log(E(", response, "))"),
                gamma = paste0("log(E(", response, "))"),
                exponential = paste0("log(E(", response, "))"),
                paste0("eta(", response, ")"))

  pieces <- character(0L)
  has_intercept <- any(terms == "(Intercept)")
  if (has_intercept) {
    pieces <- c(pieces, .memwas_format_number(beta[terms == "(Intercept)"][1L], digits = digits))
  }
  for (j in seq_along(beta)) {
    if (terms[j] == "(Intercept)") next
    bj <- beta[j]
    if (!is.finite(bj)) bj <- NA_real_
    sign_txt <- if (!is.na(bj) && bj < 0) " - " else if (length(pieces) == 0L) "" else " + "
    pieces <- c(pieces, paste0(sign_txt, .memwas_format_number(abs(bj), digits = digits), " * ", terms[j]))
  }
  rhs <- paste(pieces, collapse = "")
  if (!nzchar(rhs)) rhs <- "0"
  paste0(lhs, " = ", rhs)
}

.memwas_random_effects_table <- function(fit, settings = NULL) {
  if (!is.null(fit) && !is.null(fit$random_effects_table) &&
      is.data.frame(fit$random_effects_table)) {
    return(fit$random_effects_table)
  }
  if (is.null(fit) || is.null(fit$random_effects) || length(fit$random_effects) == 0L) {
    return(.memwas_empty_table(id = character(0L), effect = character(0L),
                               estimate = numeric(0L), std_error = numeric(0L),
                               conf_low = numeric(0L), conf_high = numeric(0L),
                               ci_level = numeric(0L)))
  }
  re <- fit$random_effects
  ids <- names(re)
  if (is.null(ids) || length(ids) != length(re)) ids <- as.character(seq_along(re))
  rows <- vector("list", length(re))
  for (i in seq_along(re)) {
    bi <- re[[i]]
    if (is.null(bi) || length(bi) == 0L) {
      rows[[i]] <- NULL
      next
    }
    eff <- names(bi)
    if (is.null(eff) || length(eff) != length(bi)) eff <- paste0("b", seq_along(bi))
    rows[[i]] <- data.frame(id = rep(ids[i], length(bi)),
                            effect = eff,
                            estimate = as.numeric(bi),
                            std_error = rep(NA_real_, length(bi)),
                            conf_low = rep(NA_real_, length(bi)),
                            conf_high = rep(NA_real_, length(bi)),
                            ci_level = rep(NA_real_, length(bi)),
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
  }
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) {
    return(.memwas_empty_table(id = character(0L), effect = character(0L),
                               estimate = numeric(0L), std_error = numeric(0L),
                               conf_low = numeric(0L), conf_high = numeric(0L),
                               ci_level = numeric(0L)))
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

.memwas_random_formula_with_coefficients <- function(fit, settings = NULL) {
  random_formula <- NULL
  id_name <- NULL
  if (!is.null(settings)) {
    random_formula <- settings$random
    id_name <- settings$id
  }
  if (is.null(random_formula) && !is.null(fit) && !is.null(fit$settings)) random_formula <- fit$settings$random
  if (is.null(id_name) && !is.null(fit) && !is.null(fit$settings)) id_name <- fit$settings$id
  rf <- .memwas_parse_random_formula(random_formula)
  rf_txt <- .memwas_deparse_one_line(rf)
  if (is.null(id_name) || length(id_name) == 0L || is.na(id_name[1L]) || !nzchar(id_name[1L])) {
    id_name <- "id"
  } else {
    id_name <- as.character(id_name[1L])
  }
  q <- if (!is.null(fit) && !is.null(fit$random_covariance)) ncol(as.matrix(fit$random_covariance)) else 0L
  n_id <- if (!is.null(fit) && !is.null(fit$random_effects)) length(fit$random_effects) else 0L
  paste0("b_", id_name, " terms: ", rf_txt, " | ", id_name,
         "; b_", id_name, " ~ N(0, D)",
         if (q > 0L) paste0(", D is ", q, " x ", q, " and stored in coefficients$random_covariance") else "",
         if (n_id > 0L) paste0("; ", n_id, " subject-specific BLUP set(s) stored in coefficients$random_effects") else "")
}

.memwas_flatten_numeric_components <- function(x, prefix = "") {
  if (is.null(x)) {
    return(.memwas_empty_table(component = character(0L), parameter = character(0L), estimate = numeric(0L)))
  }
  rows <- list()
  add_row <- function(component, parameter, estimate) {
    data.frame(component = component, parameter = parameter, estimate = as.numeric(estimate),
               stringsAsFactors = FALSE, check.names = FALSE)
  }
  walk <- function(obj, nm) {
    if (is.null(obj)) return(NULL)
    if (is.matrix(obj) && is.numeric(obj)) {
      nr <- nrow(obj)
      nc <- ncol(obj)
      rn <- rownames(obj) %||% as.character(seq_len(nr))
      cn <- colnames(obj) %||% as.character(seq_len(nc))
      rr <- vector("list", nr * nc)
      k <- 0L
      for (i in seq_len(nr)) {
        for (j in seq_len(nc)) {
          k <- k + 1L
          rr[[k]] <- add_row(nm, paste0("[", rn[i], ",", cn[j], "]"), obj[i, j])
        }
      }
      rows[[length(rows) + 1L]] <<- do.call(rbind, rr)
      return(NULL)
    }
    if (is.numeric(obj) || is.integer(obj)) {
      vals <- as.numeric(obj)
      nms <- names(obj)
      if (is.null(nms) || length(nms) != length(vals)) nms <- as.character(seq_along(vals))
      rows[[length(rows) + 1L]] <<- data.frame(component = rep(nm, length(vals)),
                                                parameter = nms,
                                                estimate = vals,
                                                stringsAsFactors = FALSE,
                                                check.names = FALSE)
      return(NULL)
    }
    if (is.list(obj)) {
      nms <- names(obj)
      if (is.null(nms) || length(nms) != length(obj)) nms <- paste0("component", seq_along(obj))
      for (i in seq_along(obj)) {
        next_nm <- if (nzchar(nm)) paste(nm, nms[i], sep = ".") else nms[i]
        walk(obj[[i]], next_nm)
      }
    }
    NULL
  }
  walk(x, prefix)
  if (!length(rows)) {
    return(.memwas_empty_table(component = character(0L), parameter = character(0L), estimate = numeric(0L)))
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

.memwas_autocorrelation_formula <- function(ac, digits = 6L) {
  if (is.null(ac) || length(ac) == 0L) return("R_i = I")
  type <- as.character(ac$type %||% "UNKNOWN")
  type_upper <- toupper(type)
  if (type_upper %in% c("NONE", "INDEPENDENCE", "INDEPENDENT")) return("R_i = I")
  if (type_upper %in% c("AR(1)", "AR1")) {
    return(paste0("Corr(e_it, e_is) = rho^|t-s|; rho = ", .memwas_format_number(ac$rho, digits)))
  }
  if (grepl("^AR\\(", type_upper)) {
    phi <- ac$ar_coefficients %||% numeric(0L)
    pacf <- ac$partial_autocorrelations %||% numeric(0L)
    return(paste0(type, "; ar_coefficients = [", .memwas_format_named_numeric(phi, digits),
                  "]; partial_autocorrelations = [", .memwas_format_named_numeric(pacf, digits), "]"))
  }
  if (type_upper %in% c("ARMA(1,1)", "ARMA11")) {
    return(paste0("ARMA(1,1); phi = ", .memwas_format_number(ac$phi, digits),
                  "; theta = ", .memwas_format_number(ac$theta, digits)))
  }
  if (type_upper %in% c("CS", "COMPOUND SYMMETRY")) {
    return(paste0("Corr(e_it, e_is) = rho for t != s; rho = ", .memwas_format_number(ac$rho, digits)))
  }
  if (type_upper %in% c("TOEP", "TOEPLITZ")) {
    return(paste0("Toeplitz residual correlation; lag_correlations = [",
                  .memwas_format_named_numeric(ac$lag_correlations %||% numeric(0L), digits), "]"))
  }
  if (type_upper %in% c("UN", "UNSTRUCTURED")) {
    return("Unstructured residual correlation; full correlation matrix stored in autocorrelation$estimates.")
  }
  paste0(type, "; parameters stored in autocorrelation$estimates")
}

.memwas_autocor_regularization_details <- function(control) {
  control <- control %||% list()
  reg <- control$autocor_regularization %||% NULL
  lambda <- control$autocor_penalty %||% 0
  type <- control$autocor_regularization_type %||% "L2"
  alpha <- control$autocor_regularization_alpha %||% 0
  enabled <- TRUE
  if (is.list(reg)) {
    enabled <- isTRUE(reg$enabled %||% TRUE)
    lambda <- reg$lambda %||% lambda
    type <- reg$type %||% type
    alpha <- reg$alpha %||% alpha
  } else if (is.numeric(reg) && length(reg) >= 1L) {
    lambda <- reg[1L]
  }
  lambda <- suppressWarnings(as.numeric(lambda)[1L])
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  if (!is.finite(lambda)) lambda <- 0
  if (!is.finite(alpha)) alpha <- 0
  alpha <- pmin(pmax(alpha, 0), 1)
  type <- tolower(as.character(type)[1L])
  if (is.na(type) || !nzchar(type)) type <- "l2"
  list(enabled = enabled, lambda = lambda, type = type, alpha = alpha)
}

.memwas_fixed_penalty_formula <- function(L1_penalty, L2_penalty, digits = 6L) {
  L1_penalty <- suppressWarnings(as.numeric(L1_penalty)[1L])
  L2_penalty <- suppressWarnings(as.numeric(L2_penalty)[1L])
  if (!is.finite(L1_penalty)) L1_penalty <- 0
  if (!is.finite(L2_penalty)) L2_penalty <- 0
  if (L1_penalty <= 0 && L2_penalty <= 0) return("P_fixed(beta) = 0")
  parts <- character(0L)
  if (L1_penalty > 0) {
    parts <- c(parts, paste0(.memwas_format_number(L1_penalty, digits),
                             " * sum_{j: beta_j is not intercept} |beta_j|"))
  }
  if (L2_penalty > 0) {
    parts <- c(parts, paste0("0.5 * ", .memwas_format_number(L2_penalty, digits),
                             " * sum_{j: beta_j is not intercept} beta_j^2"))
  }
  paste0("P_fixed(beta) = ", paste(parts, collapse = " + "))
}

.memwas_autocor_penalty_formula <- function(control, digits = 6L) {
  det <- .memwas_autocor_regularization_details(control)
  if (!isTRUE(det$enabled) || det$lambda <= 0 || det$type == "none") return("P_autocor(theta) = 0")
  lam <- .memwas_format_number(det$lambda, digits)
  alp <- .memwas_format_number(det$alpha, digits)
  if (det$type == "l1") return(paste0("P_autocor(theta) = ", lam, " * sum_k |theta_k|"))
  if (det$type == "elasticnet") {
    return(paste0("P_autocor(theta) = ", lam, " * (", alp,
                  " * sum_k |theta_k| + 0.5 * (1 - ", alp, ") * sum_k theta_k^2)"))
  }
  paste0("P_autocor(theta) = 0.5 * ", lam, " * sum_k theta_k^2")
}

.memwas_penalty_summary_for_best_model <- function(best, fit = NULL, settings = NULL, digits = 6L) {
  L1 <- if ("L1_penalty" %in% names(best) && !is.null(best$L1_penalty)) best$L1_penalty[1L] else if (!is.null(settings)) settings$L1_penalty else NA_real_
  L2 <- if ("L2_penalty" %in% names(best) && !is.null(best$L2_penalty)) best$L2_penalty[1L] else if (!is.null(settings)) settings$L2_penalty else NA_real_
  lambda <- if ("lambda" %in% names(best) && !is.null(best$lambda)) best$lambda[1L] else NA_real_
  alpha <- if ("alpha" %in% names(best) && !is.null(best$alpha)) best$alpha[1L] else NA_real_
  control <- if (!is.null(settings)) settings$control else list()
  ac_reg <- .memwas_autocor_regularization_details(control)
  fixed_value <- if (!is.null(fit) && !is.null(fit$metrics)) fit$metrics$penalty %||% NA_real_ else NA_real_
  autocor_value <- if (!is.null(fit) && !is.null(fit$metrics)) fit$metrics$autocorrelation_penalty %||% NA_real_ else NA_real_
  list(lambda = as.numeric(lambda),
       alpha = as.numeric(alpha),
       L1_penalty = as.numeric(L1),
       L2_penalty = as.numeric(L2),
       fixed_effect_formula = .memwas_fixed_penalty_formula(L1, L2, digits = digits),
       fixed_effect_value = as.numeric(fixed_value),
       autocorrelation_formula = .memwas_autocor_penalty_formula(control, digits = digits),
       autocorrelation_value = as.numeric(autocor_value),
       autocorrelation_regularization = ac_reg)
}

.memwas_build_tuned_best_model <- function(best, settings, final_fit = NULL,
                                           final_fit_error = NULL, digits = 6L) {
  best <- best[1L, , drop = FALSE]
  fit_available <- !is.null(final_fit)
  family <- if (fit_available) final_fit$family else settings$family
  method <- if (fit_available && !is.null(final_fit$metrics)) final_fit$metrics$method %||% settings$method else settings$method
  formal_formula <- if (fit_available) final_fit$formal_formula %||% final_fit$formula else settings$formal_formula %||% settings$formula
  original_formula <- settings$formula
  fixed_table <- .memwas_fixed_effect_coefficient_table(final_fit)
  random_table <- .memwas_random_effects_table(final_fit, settings)
  ac <- if (fit_available) final_fit$autocorrelation else list(type = .memwas_normalize_autocor(settings$autocor)$label)
  ac_est <- .memwas_flatten_numeric_components(ac, prefix = "autocorrelation")
  penalty <- .memwas_penalty_summary_for_best_model(best, fit = final_fit, settings = settings, digits = digits)
  dot_spec <- settings$dot_spec %||% .memwas_validate_dot_settings()

  formulas <- list(
    original = .memwas_deparse_one_line(original_formula),
    formal = .memwas_deparse_one_line(formal_formula),
    fixed_effects = .memwas_deparse_one_line(stats::delete.response(stats::terms(formal_formula))),
    fixed_effects_with_coefficients = .memwas_fixed_formula_with_coefficients(final_fit, digits = digits),
    random_effects = paste0(.memwas_deparse_one_line(.memwas_parse_random_formula(settings$random)), " | ", settings$id),
    random_effects_with_coefficients = .memwas_random_formula_with_coefficients(final_fit, settings),
    fixed_effect_penalty = penalty$fixed_effect_formula,
    autocorrelation = .memwas_autocorrelation_formula(ac, digits = digits),
    autocorrelation_penalty = penalty$autocorrelation_formula
  )

  fallback_approximation <- if (identical(family, "gaussian")) "exact_gaussian" else settings$approximation %||% "laplace"
  approximation <- if (fit_available) final_fit$approximation %||% fallback_approximation else fallback_approximation
  approximation_label <- if (fit_available) {
    final_fit$approximation_label %||% .memwas_approximation_label(approximation)
  } else {
    .memwas_approximation_label(approximation)
  }

  out <- list(
    final_fit_available = fit_available,
    final_fit_error = final_fit_error,
    family = family,
    method = method,
    approximation = approximation,
    approximation_label = approximation_label,
    init_approximation = if (fit_available) final_fit$init_approximation %||% settings$init_approximation %||% NA_character_ else settings$init_approximation %||% NA_character_,
    se_method = if (fit_available) final_fit$se_method %||% settings$se_method %||% NA_character_ else settings$se_method %||% NA_character_,
    formulas = formulas,
    coefficients = list(
      fixed_effects = fixed_table,
      random_effects = random_table,
      random_covariance = if (fit_available) final_fit$random_covariance else NULL,
      residual_sigma = if (fit_available) final_fit$residual_sigma else NA_real_
    ),
    dot_spec = dot_spec,
    directional_tests = dot_spec$table,
    penalties = penalty,
    autocorrelation = list(
      structure = ac$type %||% NA_character_,
      formula = formulas$autocorrelation,
      estimates = ac_est,
      raw = ac
    ),
    metrics = if (fit_available) final_fit$metrics else NULL,
    convergence = if (fit_available) final_fit$convergence else NA_integer_,
    selected_terms = if (nrow(fixed_table)) fixed_table$term[fixed_table$selected] else character(0L),
    n_fixed_effects = nrow(fixed_table),
    n_selected_fixed_effects = if (nrow(fixed_table)) sum(fixed_table$selected) else 0L,
    n_random_effect_rows = nrow(random_table)
  )
  class(out) <- "MEMWAS_tuned_best_model"
  out
}

.memwas_attach_best_model_columns <- function(best, best_model) {
  best <- best[1L, , drop = FALSE]
  best$final_fit_available <- isTRUE(best_model$final_fit_available)
  best$n_fixed_effects <- best_model$n_fixed_effects %||% NA_integer_
  best$n_selected_fixed_effects <- best_model$n_selected_fixed_effects %||% NA_integer_
  best$n_random_effect_rows <- best_model$n_random_effect_rows %||% NA_integer_
  best$fixed_effects_formula <- best_model$formulas$fixed_effects_with_coefficients %||% NA_character_
  best$random_effects_formula <- best_model$formulas$random_effects_with_coefficients %||% NA_character_
  best$fixed_effect_penalty_formula <- best_model$formulas$fixed_effect_penalty %||% NA_character_
  best$autocorrelation_formula <- best_model$formulas$autocorrelation %||% NA_character_
  best$autocorrelation_penalty_formula <- best_model$formulas$autocorrelation_penalty %||% NA_character_
  best
}
