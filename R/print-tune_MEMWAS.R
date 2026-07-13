#' @title Print a MEMWAS tuning result
#' @description Print the best penalty pair, selected final-model formulas, coefficient availability, and optimization diagnostics from a `tune_MEMWAS` object.
#' @param x tune_MEMWAS. Tuning object returned by `tune_MEMWAS()`.
#' @param ... list. Additional arguments; currently unused.
#' @returns invisible. Invisibly returns `x`.
#' @examples
#' \dontrun{
#' print(tuned)
#' }
#' @export
print.tune_MEMWAS <- function(x, ...) {
  cat("\n================ MEMWAS Tuning Summary ================\n")

  b <- x$best[1L, , drop = FALSE]
  bm <- x$best_model %||% NULL

  cat("\nBest hyperparameters:\n")
  cat("  lambda      :", signif(b$lambda, 6), "\n")
  cat("  log_lambda  :", signif(b$log_lambda, 6), "\n")
  cat("  alpha       :", signif(b$alpha, 6), "\n")
  cat("  L1_penalty  :", signif(b$L1_penalty, 6), "\n")
  cat("  L2_penalty  :", signif(b$L2_penalty, 6), "\n")

  cat("\nCross-validation performance:\n")
  cat("  median", x$tuning_control$metric, ":", signif(b$median_metric, 6), "\n")
  cat("  stability (", x$tuning_control$stability_metric, "):", signif(b$stability, 6), "\n", sep = "")
  cat("  selection score:", signif(b$selection_score, 6), "\n")
  cat("  successful folds:", b$n_successful_folds, "/", b$n_folds, "\n", sep = "")

  cat("\nBest final model:\n")
  if (is.null(bm)) {
    cat("  Structured best-model details are not available.\n")
  } else {
    cat("  final refit available:", if (isTRUE(bm$final_fit_available)) "yes" else "no", "\n")
    cat("  family:", bm$family %||% "NA", "\n")
    cat("  method:", bm$method %||% "NA", "\n")
    .memwas_cat_wrapped("formal formula: ", bm$formulas$formal)
    .memwas_cat_wrapped("fixed effects: ", bm$formulas$fixed_effects_with_coefficients)
    .memwas_cat_wrapped("random effects: ", bm$formulas$random_effects_with_coefficients)
    .memwas_cat_wrapped("fixed penalty: ", bm$formulas$fixed_effect_penalty)
    cat("  fixed penalty value:", .memwas_format_number(bm$penalties$fixed_effect_value), "\n")
    .memwas_cat_wrapped("autocorrelation: ", bm$formulas$autocorrelation)
    .memwas_cat_wrapped("autocor penalty: ", bm$formulas$autocorrelation_penalty)
    cat("  autocor penalty value:", .memwas_format_number(bm$penalties$autocorrelation_value), "\n")
    cat("  fixed effects stored:", bm$n_fixed_effects, "\n")
    cat("  selected/nonzero fixed effects:", bm$n_selected_fixed_effects, "\n")
    cat("  random-effect coefficient rows:", bm$n_random_effect_rows, "\n")
  }

  cat("\nOptimization:\n")
  cat("  surrogate:", x$tuning_control$surrogate, "\n")
  cat("  evaluations:", x$convergence$total_evaluations, "\n")
  cat("  completed iterations:", x$convergence$completed_iterations, "\n")
  cat("  stop reason:", x$convergence$stop_reason, "\n")

  if (!is.null(x$final_fit_error) && length(x$final_fit_error) > 0L &&
      !is.na(x$final_fit_error[1L]) && nzchar(x$final_fit_error[1L])) {
    cat("\nFinal fit note/error:\n  ", x$final_fit_error[1L], "\n", sep = "")
  }
  cat("\n======================================================\n")
  invisible(x)
}
