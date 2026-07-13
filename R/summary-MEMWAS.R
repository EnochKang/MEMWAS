#' @title Summarize a MEMWAS settings object
#' @description Create a concise summary of model setup, directional one-tailed test (DOT) settings, non-linear screening results, and GLMM assumption-screening results from `set_MEMWAS()`.
#' @param object MEMWAS. A settings object returned by `set_MEMWAS()`.
#' @param ... list. Additional arguments; currently unused.
#' @returns summary.MEMWAS. A summary object with call, family, formula, DOT settings, screening table, selected spline variables, turning points, and GLMM assumption-screening tables.
#' @examples
#' \dontrun{
#' summary(setup)
#' }
#' @export
summary.MEMWAS <- function(object, ...) {

  out <- list()

  out$call <- object$call
  out$family <- object$family
  out$approximation <- object$approximation %||% "laplace"
  out$init_approximation <- object$init_approximation %||% "variational_inference"
  out$se_method <- object$se_method %||% "hessian"
  out$approximation_label <- if (out$family == "gaussian") "Exact Gaussian marginal likelihood" else .memwas_approximation_label(out$approximation)
  out$formula <- object$formula
  out$dot_spec <- object$dot_spec %||% .memwas_validate_dot_settings(
    dot_predictors = object$dot_predictors %||% NULL,
    dot_alternative = object$dot_alternative %||% NULL,
    dot_threshold = object$dot_threshold %||% 0,
    dot_alpha = object$dot_alpha %||% 0.05
  )
  out$directional_tests <- out$dot_spec$table

  # GLMM assumption-screening tables
  out$assumption_check_spec <- .memwas_assumption_spec_from_settings(object, default = "None")
  out$assumption_check_methods <- object$assumption_check_methods %||% .memwas_assumption_spec_table(out$assumption_check_spec)
  out$assumption_checks <- object$assumption_checks %||% list(summary = .memwas_assumption_empty_table(), tables = list(), note = "No GLMM assumption-screening results available.")
  out$assumption_check_table <- .memwas_get_assumption_summary_table(out$assumption_checks)
  out$assumption_check_tables <- .memwas_get_assumption_split_tables(out$assumption_checks)
  out$assumption_check_note <- if (is.list(out$assumption_checks) && !is.null(out$assumption_checks$note)) as.character(out$assumption_checks$note)[1L] else NA_character_
  out$assumption_check_model <- object$assumption_check_model %||% if (is.list(out$assumption_checks)) out$assumption_checks$stage else NA_character_

  # Nonlinear screening table
  out$nonlinear_summary <- object$nonlinear_summary

  # Extract turning points (only for selected nonlinear vars)
  tp <- object$turning_points
  tp <- tp[object$spline_variables]

  out$turning_points <- tp
  out$spline_variables <- object$spline_variables

  class(out) <- "summary.MEMWAS"
  return(out)
}
