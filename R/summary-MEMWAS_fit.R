#' @title Summarize a fitted MEMWAS model
#' @description Summarize fixed effects with DOT alternatives and alpha levels, random-effect BLUPs with uncertainty intervals, fit metrics, GLMM assumption screening, selected non-linear effects, turning points, and segment-specific slopes from a `MEMWAS_fit` object.
#' @param object MEMWAS_fit. A fitted model object returned by `fit_MEMWAS()`.
#' @param random_effects_top_n integer. Number of subject-specific random-effect rows to preview in printed summaries.
#' @param ... list. Additional arguments; currently unused.
#' @returns summary.MEMWAS_fit. A summary object with coefficient table, DOT settings, random-effect BLUPs, metrics, GLMM assumption-screening tables, non-linear effects, and segment-specific slopes.
#' @examples
#' \dontrun{
#' summary(fit)
#' }
#' @export
summary.MEMWAS_fit <- function(object, random_effects_top_n = 10L, ...) {

  rcs_basis <- function(x, knots) {
    x <- as.numeric(x)
    knots <- sort(unique(as.numeric(knots)))
    K <- length(knots)
    if (K < 3L) stop("At least three unique knots are required for restricted cubic splines.", call. = FALSE)
    tp <- function(z, k) pmax(z - k, 0)^3
    denom <- max(knots[K] - knots[K - 1L], 1e-12)
    B <- matrix(NA_real_, nrow = length(x), ncol = K - 2L)
    for (j in seq_len(K - 2L)) {
      B[, j] <- tp(x, knots[j]) - tp(x, knots[K - 1L]) * ((knots[K] - knots[j]) / denom) +
        tp(x, knots[K]) * ((knots[K - 1L] - knots[j]) / denom)
    }
    scale <- max(diff(range(knots)), 1e-12)^3
    B <- B / scale
    colnames(B) <- paste0("spline", seq_len(ncol(B)))
    B
  }

  .compute_segment_effects <- function(var, fit) {

    if (!var %in% fit$spline_variables) return(NULL)

    info <- fit$settings$spline_info[[var]]
    knots <- info$knots
    basis_names <- info$basis_names

    co <- fit$coefficients

    beta_linear <- if (var %in% names(co)) co[var] else 0
    beta_basis <- co[basis_names]

    # numerical derivative function
    f <- function(x) {
      B <- attr(rcs_basis(x, knots), "dim") # dummy
      B <- rcs_basis(x, knots)
      beta_linear * x + as.vector(B %*% beta_basis)
    }

    # grid
    xgrid <- seq(min(knots), max(knots), length.out = 200)
    ygrid <- f(xgrid)

    dydx <- diff(ygrid) / diff(xgrid)

    segments <- cut(xgrid[-1],
                    breaks = c(-Inf, fit$settings$turning_points[[var]], Inf),
                    include.lowest = TRUE)

    seg_df <- aggregate(dydx, list(segment = segments), mean)

    names(seg_df)[2] <- "slope"

    return(seg_df)
  }

  random_effects_top_n <- as.integer(random_effects_top_n[1L])
  if (!is.finite(random_effects_top_n) || random_effects_top_n < 1L) random_effects_top_n <- 10L

  out <- list()

  out$call <- object$call
  out$formula <- object$formal_formula
  out$family <- object$family
  out$approximation <- object$approximation %||% if (object$family == "gaussian") "exact_gaussian" else "laplace"
  out$approximation_label <- object$approximation_label %||% if (object$family == "gaussian") "Exact Gaussian marginal likelihood" else .memwas_approximation_label(out$approximation)
  out$family_link <- object$family_link %||% NA_character_
  out$init_approximation <- object$init_approximation %||% object$settings$init_approximation %||% NA_character_
  out$init_approximation_label <- object$init_approximation_label %||% if (!is.na(out$init_approximation)) .memwas_approximation_label(out$init_approximation) else NA_character_
  out$se_method <- object$se_method %||% object$settings$se_method %||% NA_character_
  out$se_source <- object$se_source %||% object$metrics$se_source %||% NA_character_
  out$approximation_diagnostics <- object$approximation_diagnostics %||% list()
  out$profile_ci <- object$profile_ci %||% data.frame()
  out$family_parameters <- object$family_parameters %||% list()

  out$coefficients <- object$coefficient_table
  out$dot_spec <- object$dot_spec %||% object$settings$dot_spec %||% .memwas_validate_dot_settings()
  out$directional_tests <- out$dot_spec$table
  out$metrics <- object$metrics

  # GLMM assumption-screening tables
  settings_for_checks <- object$settings %||% list()
  out$assumption_check_spec <- object$assumption_check_spec %||% .memwas_assumption_spec_from_settings(settings_for_checks, default = "None")
  out$assumption_check_methods <- object$assumption_check_methods %||% settings_for_checks$assumption_check_methods %||% .memwas_assumption_spec_table(out$assumption_check_spec)
  stored_assumption_checks <- object$assumption_checks %||% NULL
  if (is.null(stored_assumption_checks) && !is.null(settings_for_checks$assumption_checks) &&
      length(settings_for_checks$assumption_checks) > 0L) {
    stored_assumption_checks <- settings_for_checks$assumption_checks
  }
  stored_stage <- if (is.list(stored_assumption_checks) && !is.null(stored_assumption_checks$stage)) {
    as.character(stored_assumption_checks$stage)[1L]
  } else {
    NA_character_
  }
  if (!is.null(stored_assumption_checks) && identical(stored_stage, "fitted_model")) {
    out$assumption_checks <- stored_assumption_checks
  } else if (!.memwas_assumption_is_empty_spec(out$assumption_check_spec)) {
    assumption_checks_try <- try(
      .memwas_run_assumption_checks(
        object,
        settings = settings_for_checks,
        spec = out$assumption_check_spec,
        stage = "fitted_model"
      ),
      silent = TRUE
    )
    out$assumption_checks <- if (inherits(assumption_checks_try, "try-error")) {
      .memwas_assumption_checks_unavailable(
        out$assumption_check_spec,
        settings = settings_for_checks,
        reason = paste("GLMM assumption screening failed:", as.character(assumption_checks_try)),
        stage = "fitted_model"
      )
    } else {
      assumption_checks_try
    }
  } else {
    out$assumption_checks <- stored_assumption_checks %||% list(
      stage = "fitted_model",
      alpha = .memwas_assumption_alpha(settings_for_checks),
      methods = out$assumption_check_spec,
      summary = .memwas_assumption_empty_table(),
      tables = list(),
      note = "No GLMM assumption screening methods were requested."
    )
  }
  out$assumption_check_table <- .memwas_get_assumption_summary_table(out$assumption_checks)
  out$assumption_check_tables <- .memwas_get_assumption_split_tables(out$assumption_checks)
  out$assumption_check_note <- if (is.list(out$assumption_checks) && !is.null(out$assumption_checks$note)) as.character(out$assumption_checks$note)[1L] else NA_character_
  out$assumption_check_model <- if (is.list(out$assumption_checks) && !is.null(out$assumption_checks$stage)) as.character(out$assumption_checks$stage)[1L] else "fitted_model"

  out$random_effects <- object$random_effects_table %||% data.frame()
  out$random_effects_preview <- .memwas_first_rows(out$random_effects, random_effects_top_n)
  out$random_effects_n_rows <- if (is.data.frame(out$random_effects)) nrow(out$random_effects) else 0L
  out$random_effects_ci_level <- object$random_effects_ci_level %||% NA_real_
  out$random_effects_note <- object$random_effects_note %||% NA_character_
  out$random_covariance <- object$random_covariance
  out$residual_sigma <- object$residual_sigma
  out$random_effects_top_n <- random_effects_top_n

  out$turning_points <- object$settings$turning_points
  out$spline_variables <- object$spline_variables

  # Segment-specific effects
  seg_list <- list()

  for (v in object$spline_variables) {
    seg_list[[v]] <- .compute_segment_effects(v, object)
  }

  out$segment_effects <- seg_list

  class(out) <- "summary.MEMWAS_fit"
  return(out)
}
