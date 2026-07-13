#' Diagnose MEMWAS non-Gaussian approximation quality
#'
#' @description Extract approximation diagnostics from a `MEMWAS_fit` or
#' `tune_MEMWAS` object. The diagnostic table reports the final approximation,
#' initialization method, random-effect dimension, convergence, Hessian status,
#' quadrature fallback, standard-error source, and subject-level approximation
#' details when available.
#' @param object A `MEMWAS_fit` object, a `tune_MEMWAS` object, or an object
#' containing `best_fit`.
#' @param ... Reserved for future extensions.
#' @returns A list of class `MEMWAS_approximation_diagnostics` containing a
#' one-row summary and optional subject-level details.
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_MEMWAS(setup)
#' diagnose_approximation(fit)
#' }
diagnose_approximation <- function(object, ...) {
  fit <- object
  if (inherits(object, "tune_MEMWAS")) fit <- object$best_fit
  if (is.null(fit) || !inherits(fit, "MEMWAS_fit")) {
    stop("`object` must be a MEMWAS_fit object or a tune_MEMWAS object with a fitted `best_fit`.", call. = FALSE)
  }

  d <- fit$approximation_diagnostics %||% list()
  details <- d$details %||% fit$approximation_details$details %||% data.frame()
  n_groups <- length(fit$groups %||% list())
  q <- ncol(fit$Z %||% matrix(numeric(0L), 0L, 0L))
  summary <- data.frame(
    family = fit$family %||% NA_character_,
    link = fit$family_link %||% NA_character_,
    approximation = fit$approximation %||% NA_character_,
    approximation_label = fit$approximation_label %||% NA_character_,
    init_approximation = fit$init_approximation %||% NA_character_,
    init_approximation_label = fit$init_approximation_label %||% NA_character_,
    nobs = fit$metrics$nobs %||% length(fit$y %||% numeric(0L)),
    n_groups = n_groups,
    random_effect_dimension = q,
    converged = isTRUE((fit$convergence %||% 1L) == 0L),
    optimizer_convergence = fit$convergence %||% NA_integer_,
    hessian_positive_definite = isTRUE(fit$metrics$hessian_positive_definite %||% d$hessian_positive_definite %||% FALSE),
    se_source = fit$se_source %||% fit$metrics$se_source %||% NA_character_,
    fallback_groups = fit$metrics$approximation_fallback_groups %||% d$fallback_groups %||% NA_integer_,
    quadrature_nodes = fit$approximation_details$n_quadrature_nodes %||% d$quadrature_nodes %||% NA_integer_,
    logLik = fit$metrics$approximate_marginal_logLik %||% fit$metrics$logLik %||% NA_real_,
    AIC = fit$metrics$AIC_approximation %||% fit$metrics$AIC %||% NA_real_,
    BIC = fit$metrics$BIC_approximation %||% fit$metrics$BIC %||% NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- list(summary = summary,
              details = details,
              random_effect_dimension = q,
              note = fit$approximation_details$note %||% fit$methodological_note %||% "")
  class(out) <- "MEMWAS_approximation_diagnostics"
  out
}

#' @export
print.MEMWAS_approximation_diagnostics <- function(x, ...) {
  cat("MEMWAS approximation diagnostics\n")
  print(x$summary, row.names = FALSE)
  if (is.data.frame(x$details) && nrow(x$details) > 0L) {
    cat("\nSubject-level approximation details: ", nrow(x$details), " group(s).\n", sep = "")
    print(utils::head(x$details, 10L), row.names = FALSE)
  }
  if (nzchar(x$note %||% "")) cat("\nNote: ", x$note, "\n", sep = "")
  invisible(x)
}

#' Compare MEMWAS approximation methods
#'
#' @description Refit the same MEMWAS model with several non-Gaussian
#' approximation strategies and return coefficient-level and fit-level
#' sensitivity summaries. This is useful for checking whether scientific
#' conclusions, directional one-tailed test decisions, or standard errors are
#' sensitive to the final approximation.
#' @param object A `MEMWAS` settings object, a `MEMWAS_fit` object, a formula,
#' or `NULL` when formula/data/id/time/random are supplied through `...`.
#' @param approximations Character vector of final approximation methods to compare.
#' @param init_approximation Initial approximation used for all refits. Defaults
#' to `"variational_inference"`.
#' @param se_method Standard-error method used in the refits. Defaults to `"hessian"`.
#' @param refit logical. Whether to refit models. Currently must be `TRUE`.
#' @param ... Arguments passed to `fit_MEMWAS()` when `object` is a formula or
#' `NULL`, or override arguments passed to `fit_MEMWAS()` for each MEMWAS
#' settings refit.
#' @returns A list of class `MEMWAS_approximation_comparison` with refitted
#' models, fit-level metrics, coefficient tables in long form, and per-term
#' sensitivity ranges.
#' @export
#' @examples
#' \dontrun{
#' cmp <- compare_approximations(
#'   setup,
#'   approximations = c("laplace", "adaptive_gauss_hermite_quadrature")
#' )
#' print(cmp)
#' }
compare_approximations <- function(object = NULL,
                                   approximations = c("laplace", "adaptive_gauss_hermite_quadrature",
                                                      "saddlepoint", "skew_corrected_laplace"),
                                   init_approximation = "variational_inference",
                                   se_method = "hessian",
                                   refit = TRUE,
                                   ...) {
  if (!isTRUE(refit)) stop("`refit = FALSE` is not currently supported.", call. = FALSE)
  approximations <- unique(unname(vapply(approximations, .memwas_normalize_approximation, character(1L))))
  init_approximation <- .memwas_normalize_init_approximation(init_approximation)
  se_method <- .memwas_normalize_se_method(se_method)
  dots <- list(...)

  base_settings <- NULL
  if (inherits(object, "MEMWAS_fit")) base_settings <- object$settings
  if (inherits(object, "MEMWAS")) base_settings <- object

  fits <- vector("list", length(approximations))
  names(fits) <- approximations
  errors <- setNames(rep(NA_character_, length(approximations)), approximations)

  for (a in approximations) {
    if (!is.null(base_settings)) {
      settings <- base_settings
      settings$approximation <- a
      settings$init_approximation <- init_approximation
      settings$se_method <- se_method
      fit_call <- c(list(object = settings, verbose = FALSE), dots)
    } else {
      fit_call <- c(list(object = object, approximation = a,
                         init_approximation = init_approximation,
                         se_method = se_method, verbose = FALSE), dots)
    }
    fit_i <- try(do.call(fit_MEMWAS, fit_call), silent = TRUE)
    if (inherits(fit_i, "try-error")) {
      errors[a] <- as.character(fit_i)
    } else {
      fits[[a]] <- fit_i
    }
  }

  fit_rows <- lapply(names(fits), function(a) {
    f <- fits[[a]]
    if (is.null(f)) {
      return(data.frame(approximation = a, converged = FALSE, logLik = NA_real_,
                        AIC = NA_real_, BIC = NA_real_, se_source = NA_character_,
                        fallback_groups = NA_integer_, error = errors[[a]],
                        stringsAsFactors = FALSE, check.names = FALSE))
    }
    data.frame(approximation = a,
               approximation_label = f$approximation_label %||% .memwas_approximation_label(a),
               init_approximation = f$init_approximation %||% init_approximation,
               converged = isTRUE((f$convergence %||% 1L) == 0L),
               logLik = f$metrics$approximate_marginal_logLik %||% f$metrics$logLik %||% NA_real_,
               AIC = f$metrics$AIC_approximation %||% f$metrics$AIC %||% NA_real_,
               BIC = f$metrics$BIC_approximation %||% f$metrics$BIC %||% NA_real_,
               se_source = f$se_source %||% f$metrics$se_source %||% NA_character_,
               fallback_groups = f$metrics$approximation_fallback_groups %||% NA_integer_,
               error = NA_character_, stringsAsFactors = FALSE, check.names = FALSE)
  })
  fit_metrics <- do.call(rbind, fit_rows)
  row.names(fit_metrics) <- NULL

  coef_rows <- lapply(names(fits), function(a) {
    f <- fits[[a]]
    tab <- if (!is.null(f)) f$coefficient_table else NULL
    if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L) return(NULL)
    tab$approximation <- a
    keep <- intersect(c("approximation", "term", "estimate", "std_error", "statistic",
                        "p_value", "alternative", "alpha", "significant"), names(tab))
    tab[, keep, drop = FALSE]
  })
  coefficients_long <- do.call(rbind, coef_rows)
  if (is.null(coefficients_long)) coefficients_long <- data.frame()
  row.names(coefficients_long) <- NULL

  sensitivity <- data.frame()
  if (is.data.frame(coefficients_long) && nrow(coefficients_long) > 0L && "estimate" %in% names(coefficients_long)) {
    split_terms <- split(coefficients_long, coefficients_long$term)
    sensitivity <- do.call(rbind, lapply(split_terms, function(z) {
      est <- z$estimate[is.finite(z$estimate)]
      se <- if ("std_error" %in% names(z)) z$std_error[is.finite(z$std_error)] else numeric(0L)
      data.frame(term = z$term[1L],
                 n_approximations = length(unique(z$approximation)),
                 min_estimate = if (length(est)) min(est) else NA_real_,
                 max_estimate = if (length(est)) max(est) else NA_real_,
                 max_abs_estimate_difference = if (length(est)) max(est) - min(est) else NA_real_,
                 min_std_error = if (length(se)) min(se) else NA_real_,
                 max_std_error = if (length(se)) max(se) else NA_real_,
                 stringsAsFactors = FALSE, check.names = FALSE)
    }))
    row.names(sensitivity) <- NULL
  }

  out <- list(approximations = approximations,
              init_approximation = init_approximation,
              se_method = se_method,
              fits = fits,
              errors = errors,
              fit_metrics = fit_metrics,
              coefficients = coefficients_long,
              sensitivity = sensitivity)
  class(out) <- "MEMWAS_approximation_comparison"
  out
}

#' @export
print.MEMWAS_approximation_comparison <- function(x, ...) {
  cat("MEMWAS approximation sensitivity comparison\n")
  cat("Initial approximation: ", .memwas_approximation_label(x$init_approximation), "\n", sep = "")
  cat("Standard-error method: ", paste(x$se_method, collapse = ", "), "\n", sep = "")
  cat("\nFit-level metrics:\n")
  print(x$fit_metrics, row.names = FALSE)
  if (is.data.frame(x$sensitivity) && nrow(x$sensitivity) > 0L) {
    cat("\nCoefficient sensitivity by term:\n")
    print(x$sensitivity, row.names = FALSE)
  }
  failed <- names(x$errors)[!is.na(x$errors) & nzchar(x$errors)]
  if (length(failed)) {
    cat("\nFailed approximation(s): ", paste(failed, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
