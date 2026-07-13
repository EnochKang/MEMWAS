#' @title Print a MEMWAS setup summary
#' @description Print the summary produced by `summary.MEMWAS()`.
#' @param x summary.MEMWAS. Summary object to print.
#' @param ... list. Additional arguments; currently unused.
#' @returns invisible. Invisibly returns `x`.
#' @examples
#' \dontrun{
#' print(summary(setup))
#' }
#' @export
print.summary.MEMWAS <- function(x, ...) {

  cat("\n================ MEMWAS Setup Summary ================\n")

  cat("\nCall:\n")
  print(x$call)

  cat("\nModel:\n")
  cat("  Family :", x$family, "\n")
  cat("  Final approximation:", x$approximation_label %||% x$approximation %||% "NA", "\n")
  if (!is.null(x$init_approximation)) cat("  Initialization:", .memwas_approximation_label(x$init_approximation), "\n")
  if (!is.null(x$se_method)) cat("  SE method:", paste(x$se_method, collapse = ", "), "\n")
  cat("  Formula:", deparse(x$formula), "\n")

  cat("\n--- Directional One-tailed Test Settings ---\n")

  if (!is.null(x$directional_tests) && is.data.frame(x$directional_tests) && nrow(x$directional_tests) > 0) {
    print(x$directional_tests, row.names = FALSE)
  } else if (!is.null(x$dot_spec)) {
    cat("Default for fitted coefficient terms: alternative =", x$dot_spec$default_alternative,
        ", dot_threshold =", x$dot_spec$default_dot_threshold %||% x$dot_spec$default_threshold,
        ", minimal important value =", x$dot_spec$default_minimal_important_value %||% abs(x$dot_spec$default_threshold),
        ", test threshold =", x$dot_spec$default_test_threshold %||% x$dot_spec$default_threshold,
        ", alpha =", x$dot_spec$default_alpha, "\n")
  } else {
    cat("Default two-sided Wald tests with dot_threshold 0 and alpha 0.05.\n")
  }

  cat("\n--- GLMM Assumption Screening Methods ---\n")
  if (!is.null(x$assumption_check_methods) && is.data.frame(x$assumption_check_methods) &&
      nrow(x$assumption_check_methods) > 0L) {
    print(x$assumption_check_methods, row.names = FALSE)
  } else {
    cat("No GLMM assumption-screening methods stored.\n")
  }

  cat("\n--- GLMM Assumption Screening Results ---\n")
  if (!is.null(x$assumption_check_model) && !is.na(x$assumption_check_model)) {
    cat("Screening model:", x$assumption_check_model, "\n")
  }
  .memwas_print_assumption_results(x$assumption_check_table)
  if (!is.null(x$assumption_check_note) && !is.na(x$assumption_check_note) && nzchar(x$assumption_check_note)) {
    cat("Note:", x$assumption_check_note, "\n")
  }

  cat("\n--- Non-linear Screening Results ---\n")

  if (!is.null(x$nonlinear_summary) && nrow(x$nonlinear_summary) > 0) {
    print(x$nonlinear_summary, row.names = FALSE)
  } else {
    cat("No screening results available.\n")
  }

  cat("\n--- Selected Non-linear Predictors ---\n")

  if (length(x$spline_variables) == 0) {
    cat("None\n")
  } else {
    for (v in x$spline_variables) {
      cat("\nVariable:", v, "\n")
      cat("Turning points:", paste(signif(x$turning_points[[v]], 5), collapse = ", "), "\n")
    }
  }

  cat("\n=====================================================\n")
}


# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# S3 METHOD: print summary for MEMWAS_fit objects from function MEMWAS  -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
