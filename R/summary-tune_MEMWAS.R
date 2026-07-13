#' @title Summarize MEMWAS tuning results
#' @description Summarize the best penalty pair, top observed results, final best-model formulas with coefficients and DOT columns, random-effect BLUPs with uncertainty intervals, tuning controls, and convergence diagnostics.
#' @param object tune_MEMWAS. Tuning object returned by `tune_MEMWAS()`.
#' @param top_n integer. Number of top observed hyperparameter rows to include.
#' @param random_effects_top_n integer. Number of subject-specific random-effect rows to preview.
#' @param ... list. Additional arguments; currently unused.
#' @returns summary.tune_MEMWAS. Summary object containing the best row, top results, final best-model formulas and coefficients, DOT settings, full random-effect BLUP table, controls, convergence, and final-fit status.
#' @examples
#' \dontrun{
#' summary(tuned, top_n = 5)
#' }
#' @export
summary.tune_MEMWAS <- function(object, top_n = 10L, random_effects_top_n = 10L, ...) {
  top_n <- as.integer(top_n[1L])
  if (!is.finite(top_n) || top_n < 1L) top_n <- 10L
  random_effects_top_n <- as.integer(random_effects_top_n[1L])
  if (!is.finite(random_effects_top_n) || random_effects_top_n < 1L) random_effects_top_n <- 10L

  n_results <- if (is.null(object$results)) 0L else nrow(object$results)
  top_results <- if (n_results > 0L) object$results[seq_len(min(n_results, top_n)), , drop = FALSE] else data.frame()
  bm <- object$best_model %||% NULL

  fixed_effects <- if (!is.null(bm)) bm$coefficients$fixed_effects else data.frame()
  random_effects <- if (!is.null(bm)) bm$coefficients$random_effects else data.frame()
  random_preview <- .memwas_first_rows(random_effects, random_effects_top_n)
  dot_spec <- if (!is.null(bm)) {
    bm$dot_spec %||% object$best_settings$dot_spec %||% object$settings$dot_spec %||% .memwas_validate_dot_settings()
  } else {
    object$best_settings$dot_spec %||% object$settings$dot_spec %||% .memwas_validate_dot_settings()
  }

  fallback_approximation <- if (!is.null(object$settings) && identical(object$settings$family, "gaussian")) {
    "exact_gaussian"
  } else {
    object$settings$approximation %||% "laplace"
  }
  approximation <- if (!is.null(bm)) bm$approximation %||% fallback_approximation else fallback_approximation
  approximation_label <- if (!is.null(bm)) {
    bm$approximation_label %||% .memwas_approximation_label(approximation)
  } else {
    .memwas_approximation_label(approximation)
  }

  out <- list(
    call = object$call,
    best = object$best,
    top_results = top_results,
    best_model = bm,
    approximation = approximation,
    approximation_label = approximation_label,
    init_approximation = if (!is.null(bm)) bm$init_approximation %||% object$settings$init_approximation %||% NA_character_ else object$settings$init_approximation %||% NA_character_,
    se_method = object$settings$se_method %||% if (!is.null(object$best_fit)) object$best_fit$se_method else NA_character_,
    formulas = if (!is.null(bm)) bm$formulas else NULL,
    fixed_effects = fixed_effects,
    dot_spec = dot_spec,
    directional_tests = dot_spec$table,
    random_effects = random_effects,
    random_effects_preview = random_preview,
    random_effects_n_rows = if (is.data.frame(random_effects)) nrow(random_effects) else 0L,
    random_covariance = if (!is.null(bm)) bm$coefficients$random_covariance else NULL,
    residual_sigma = if (!is.null(bm)) bm$coefficients$residual_sigma else NA_real_,
    penalties = if (!is.null(bm)) bm$penalties else NULL,
    autocorrelation = if (!is.null(bm)) bm$autocorrelation else NULL,
    metrics = if (!is.null(bm)) bm$metrics else NULL,
    tuning_control = object$tuning_control,
    convergence = object$convergence,
    final_fit_available = !is.null(object$best_fit),
    final_fit_error = object$final_fit_error,
    top_n = top_n,
    random_effects_top_n = random_effects_top_n
  )
  class(out) <- "summary.tune_MEMWAS"
  out
}
