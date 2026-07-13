#' @title Fit a MEMWAS mixed-effects model
#' @description Fit the MEMWAS mixed-effects model specified by a `MEMWAS` settings object or by direct formula/data arguments, with optional fixed-effect elastic-net penalties, random effects, residual autocorrelation, and stored spline terms.
#' @param object MEMWAS or formula or NULL. A settings object returned by `set_MEMWAS()`, a model formula, or `NULL` when all required settings are supplied manually.
#' @param formula formula. Fixed-effect model formula; ignored when `object` is a populated `MEMWAS` settings object unless explicitly supplied as an override.
#' @param family character or family. Response family. Supported options are `"gaussian"`, `"binomial"`, `"poisson"`, `"negative_binomial"`, `"gamma"`, and `"exponential"` plus aliases.
#' @param data data.frame. Long-format analysis data.
#' @param id character. Name of the subject identifier column.
#' @param time character. Name of the measurement-time column.
#' @param random formula or character. Random-effect design, such as `~ 1` or `~ 1 + time`.
#' @param autocor character. One legacy residual autocorrelation structure: `"NONE"`, `"AR(1)"`, `"AR(p)"`, `"ARMA(1,1)"`, `"CS"`, `"TOEP"`, or `"UN"`. A fully named character vector may also define multiple components.
#' @param serial Optional serial-process component or named list/vector of components created with `serial_component()`. Each component contributes `diag(h_k) u_k` and has its own covariance structure.
#' @param L1_penalty numeric. Non-negative L1 fixed-effect penalty; the intercept is unpenalized.
#' @param L2_penalty numeric. Non-negative L2 fixed-effect penalty; the intercept is unpenalized.
#' @param control list. Fitting controls, including optimizer controls, approximation controls, coordinate-descent controls, spline grid size, and autocorrelation regularization settings. Family-specific entries include `negative_binomial_theta`, `gamma_shape`, `aghq_nodes`, `aghq_max_dim`, `bootstrap_B`, and profile-likelihood controls.
#' @param method character. `"ML"` or `"REML"` for Gaussian models; non-Gaussian models use `"ML"` with the selected approximation.
#' @param random_cov character. Random-effect covariance structure: `"diagonal"` or `"unstructured"`.
#' @param approximation character. Final non-Gaussian approximation strategy. The default is `"laplace"` after variational-inference initialization. Supported options are `"laplace"`, `"variational_inference"`, `"adaptive_gauss_hermite_quadrature"`, `"adaptive_gaussian_quadrature"`, `"saddlepoint"`, `"skew_corrected_laplace"`, and `"pql"`.
#' @param init_approximation character. Initial approximation used to generate starting values for non-Gaussian final optimization. Defaults to `"variational_inference"`.
#' @param se_method character. Standard-error method for non-Gaussian final fits. Use `"hessian"`, `"cluster_sandwich"`, `"parametric_bootstrap"`, `"profile"`, or a character vector combining these options.
#' @param dot_predictors character. Optional fixed-effect coefficient terms for directional one-tailed tests (DOT). Names should match rows in the fixed-effect coefficient table, such as `"x1"` or `"groupB"`.
#' @param dot_alternative character, named character, list, or data.frame. DOT alternative(s): `"greater"` tests `beta > abs(dot_threshold)`, `"less"` tests `beta < -abs(dot_threshold)`, and `"two.sided"` uses the usual two-sided Wald test when `dot_threshold = 0` or a two-sided minimum-effect test of `|beta| > abs(dot_threshold)` when `dot_threshold > 0`. Named values override specific terms; a scalar value applies to `dot_predictors`, or to all non-intercept fixed-effect coefficient terms when `dot_predictors` is `NULL`. Data frames may be supplied to `dot_alternative`, `dot_threshold`, or `dot_alpha` with a `term`/`predictor` column plus relevant `alternative`, `threshold`/`dot_threshold`, or `alpha` columns.
#' @param dot_threshold numeric or data.frame. Predictor-specific minimally important value(s). Values are interpreted as magnitudes and converted with `abs()`: `"greater"` tests against `+abs(dot_threshold)`, `"less"` tests against `-abs(dot_threshold)`, and `"two.sided"` with a nonzero value tests whether the absolute effect exceeds that magnitude. Use a scalar, a vector aligned with `dot_predictors`, a named vector/list keyed by coefficient or model term, or a data frame with `term`/`predictor` and `threshold`/`dot_threshold` columns. Defaults to 0.
#' @param dot_alpha numeric. Significance level(s) for DOT. Use a scalar, a vector aligned with `dot_predictors`, or a named vector/list keyed by coefficient term. Defaults to 0.05.
#' @param spline_variables character. Names of predictors whose stored spline bases should be included.
#' @param spline_info list. Stored spline knot and basis-name metadata, usually produced by `set_MEMWAS()`.
#' @param turning_points list. Optional turning-point metadata from spline screening.
#' @param nonlinear_summary data.frame. Optional spline-screening summary table.
#' @param all_screened_spline_info list. Optional metadata for all screened spline predictors.
#' @param baseline_screen_metrics list. Optional fit metrics from the baseline screening model.
#' @param engine character. Execution engine: `"R"` uses base-R numerical helpers; `"cpp"` uses native C++ helpers for available numerical operations including safe exponentials, diagonals, linear solves, Cholesky factorizations, penalty kernels, coordinate descent, lag matrices, AR helper calculations, spline bases, and other calculation kernels.
#' @param verbose logical. Whether to print progress messages.
#' @param ... list. Additional overrides or metadata.
#' @details
#' `fit_MEMWAS()` first merges a settings object with any explicit user overrides and validates formula, data, subject, time, penalty, method, random covariance, and engine options. It then constructs fixed- and random-effect design matrices from R formulas and fits a marginal mixed-effects model.
#'
#' For `family = "gaussian"`, the function optimizes the marginal Gaussian likelihood, or REML when `method = "REML"` and no fixed-effect penalties are used. For `family = "binomial"`, `"poisson"`, `"negative_binomial"`, `"gamma"`, and `"exponential"`, the function uses variational inference as the default initializer and then evaluates the requested final non-Gaussian marginal approximation with native C++ subject-level kernels. Negative-binomial and gamma dispersion parameters are treated as likelihood parameters unless disabled in `control`. Fixed effects can be estimated by ridge or elastic-net style penalties via `L1_penalty` and `L2_penalty`; if an L1 penalty is present, coordinate descent is used for the fixed-effect subproblem.
#'
#' Residual correlation is selected with `autocor`. `"AR(1)"`, higher-order AR, ARMA(1,1), compound-symmetry, Toeplitz, and unstructured structures are available. Random-effect covariance is controlled by `random_cov`. The `engine` argument changes only the calculation backend: formula parsing, S3 output, and model assembly stay in R, while numerical helper routines are routed to C++ when `engine = "cpp"`.
#' @returns MEMWAS_fit. A fitted model object with coefficients, standard errors, fitted values, residuals, random effects, covariance estimates, autocorrelation estimates, fit metrics, convergence information, base-R GLMM assumption-screening tables, and original settings.
#' @examples
#' \dontrun{
#' dat <- simulate_panel_data(n_id = 100, n_time = 5)
#' fit <- fit_MEMWAS(var_1 ~ var_2 + var_3, family = "gaussian", data = dat,
#'                   id = "id", time = "time", random = ~ 1,
#'                   autocor = "AR(1)", engine = "R")
#' summary(fit)
#' }
#' @export
fit_MEMWAS <- function(
    object = NULL,
    formula = NULL,
    family = NULL,
    data = NULL, id = NULL, time = NULL,
    random = NULL,
    autocor = NULL,
    L1_penalty = NULL, L2_penalty = NULL,
    control = NULL,
    method = NULL,
    random_cov = NULL,
    approximation = NULL,
    init_approximation = NULL,
    se_method = NULL,
    dot_predictors = NULL, dot_alternative = NULL,
    dot_threshold = 0, dot_alpha = 0.05,
    spline_variables = NULL, spline_info = NULL,
    turning_points = NULL, nonlinear_summary = NULL,
    all_screened_spline_info = NULL,
    baseline_screen_metrics = NULL,
    engine = NULL,
    verbose = NULL, ...) {
  dots <- list(...)
  user_call <- match.call()

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 1: Build and validate the MEMWAS settings object, allowing manual overrides -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  supplied <- list(
    formula = !missing(formula),
    family = !missing(family),
    data = !missing(data),
    id = !missing(id),
    time = !missing(time),
    random = !missing(random),
    autocor = !missing(autocor),
    L1_penalty = !missing(L1_penalty),
    L2_penalty = !missing(L2_penalty),
    control = !missing(control),
    method = !missing(method),
    random_cov = !missing(random_cov),
    approximation = !missing(approximation),
    init_approximation = !missing(init_approximation),
    se_method = !missing(se_method),
    dot_predictors = !missing(dot_predictors),
    dot_alternative = !missing(dot_alternative),
    dot_threshold = !missing(dot_threshold),
    dot_alpha = !missing(dot_alpha),
    spline_variables = !missing(spline_variables),
    spline_info = !missing(spline_info),
    turning_points = !missing(turning_points),
    nonlinear_summary = !missing(nonlinear_summary),
    all_screened_spline_info = !missing(all_screened_spline_info),
    baseline_screen_metrics = !missing(baseline_screen_metrics),
    engine = !missing(engine),
    verbose = !missing(verbose)
  )

  message("01. Checking object from function set_MEMWAS...")

  message("02. Setting parameters...")

  if (missing(object) || is.null(object)) {
    message("MEMWAS: MEMWAS class object is found, and the model will be carried out according to the arguments for the function fit_MEMWAS.")
    settings <- .memwas_make_default_settings(user_call, verbose = TRUE)
  } else if (inherits(object, "MEMWAS")) {
    settings <- object
  } else if (inherits(object, "formula")) {
    if (supplied$formula && !is.null(formula)) {
      stop("MEMWAS: Supply the model formula either as `object` or `formula`, not both.", call. = FALSE)
    }
    settings <- .memwas_make_default_settings(user_call, verbose = TRUE)
    settings$formula <- object
    supplied$formula <- FALSE
  } else {
    stop("MEMWAS: `object` must be a MEMWAS class object produced by `set_MEMWAS()`, a model formula, or NULL when all required settings are supplied manually.",
         call. = FALSE)
  }

  if (is.null(settings$extra) || !is.list(settings$extra)) settings$extra <- list()
  manual_overrides <- names(supplied)[vapply(supplied, isTRUE, logical(1L))]

  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("", length(dots))
    override_from_dots <- which(nzchar(dot_names) & dot_names %in% names(settings))
    if (length(override_from_dots) > 0L) {
      for (i in override_from_dots) settings[[dot_names[i]]] <- dots[[i]]
      manual_overrides <- unique(c(manual_overrides, dot_names[override_from_dots]))
    }
    extra_from_dots <- setdiff(seq_along(dots), override_from_dots)
    if (length(extra_from_dots) > 0L) settings$extra <- c(settings$extra, dots[extra_from_dots])
  }

  if (supplied$formula) settings$formula <- formula
  if (supplied$family) settings$family <- family
  if (supplied$data) settings$data <- data
  if (supplied$id) settings$id <- id
  if (supplied$time) settings$time <- time
  if (supplied$random) settings$random <- random
  if (supplied$autocor) settings$autocor <- autocor
  if (supplied$L1_penalty) settings$L1_penalty <- L1_penalty
  if (supplied$L2_penalty) settings$L2_penalty <- L2_penalty
  if (supplied$control) settings$control <- control
  if (supplied$method) settings$method <- method
  if (supplied$random_cov) settings$random_cov <- random_cov
  if (supplied$approximation) settings$approximation <- approximation
  if (supplied$init_approximation) settings$init_approximation <- init_approximation
  if (supplied$se_method) settings$se_method <- se_method
  if (supplied$dot_predictors) settings$dot_predictors <- dot_predictors
  if (supplied$dot_alternative) settings$dot_alternative <- dot_alternative
  if (supplied$dot_threshold) settings$dot_threshold <- dot_threshold
  if (supplied$dot_alpha) settings$dot_alpha <- dot_alpha
  if (supplied$spline_variables) settings$spline_variables <- spline_variables
  if (supplied$spline_info) settings$spline_info <- spline_info
  if (supplied$turning_points) settings$turning_points <- turning_points
  if (supplied$nonlinear_summary) settings$nonlinear_summary <- nonlinear_summary
  if (supplied$all_screened_spline_info) settings$all_screened_spline_info <- all_screened_spline_info
  if (supplied$baseline_screen_metrics) settings$baseline_screen_metrics <- baseline_screen_metrics
  if (supplied$engine) settings$engine <- engine
  if (supplied$verbose) settings$verbose <- verbose

  required <- c("formula", "family", "data", "id", "time", "L1_penalty", "L2_penalty")
  missing_required <- required[vapply(required, function(z) !z %in% names(settings) || is.null(settings[[z]]), logical(1))]
  if (length(missing_required) > 0L) {
    stop("MEMWAS: MEMWAS settings are missing required setting(s): ",
         paste(missing_required, collapse = ", "),
         ". Supply a `set_MEMWAS()` object or pass these arguments directly to `fit_MEMWAS()`.",
         call. = FALSE)
  }
  if (is.null(settings$random)) settings$random <- stats::as.formula("~ 1")
  if (is.null(settings$method)) settings$method <- "ML"
  if (is.null(settings$random_cov)) settings$random_cov <- "diagonal"
  if (is.null(settings$approximation)) settings$approximation <- "laplace"
  if (is.null(settings$init_approximation)) settings$init_approximation <- "variational_inference"
  if (is.null(settings$se_method)) settings$se_method <- "hessian"
  if (is.null(settings$control)) settings$control <- list()
  if (is.null(settings$dot_threshold)) settings$dot_threshold <- 0
  if (is.null(settings$dot_alpha)) settings$dot_alpha <- 0.05
  if (is.null(settings$spline_variables)) settings$spline_variables <- character(0L)
  if (is.null(settings$spline_info)) settings$spline_info <- list()
  if (is.null(settings$turning_points)) settings$turning_points <- list()
  if (is.null(settings$nonlinear_summary)) settings$nonlinear_summary <- data.frame()
  if (is.null(settings$all_screened_spline_info)) settings$all_screened_spline_info <- list()
  if (is.null(settings$engine)) settings$engine <- "R"
  settings$engine <- .memwas_normalize_engine(settings$engine)
  if (is.null(settings$verbose)) settings$verbose <- TRUE

  if (!inherits(settings$formula, "formula")) stop("`formula` must be an R formula.", call. = FALSE)
  if (!is.data.frame(settings$data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(settings$id) || length(settings$id) != 1L || !settings$id %in% names(settings$data)) {
    stop("MEMWAS: `id` must be the name of an existing subject-id column in `data`.", call. = FALSE)
  }
  if (!is.character(settings$time) || length(settings$time) != 1L || !settings$time %in% names(settings$data)) {
    stop("MEMWAS: `time` must be the name of an existing measurement-time column in `data`.", call. = FALSE)
  }
  if (!is.list(settings$control)) stop("MEMWAS: `control` must be a list.", call. = FALSE)
  if (!is.list(settings$spline_info)) stop("MEMWAS: `spline_info` must be a list.", call. = FALSE)

  settings$manual_overrides <- manual_overrides
  if (is.null(settings$call)) settings$call <- user_call
  if (!inherits(settings, "MEMWAS")) class(settings) <- c("MEMWAS", class(settings))


  backend <- settings$engine
  verbose <- isTRUE(settings$verbose)

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 2: Use internal helper functions for base-R and C++-backed estimation -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 3: Normalize options and prepare spline-augmented formal formula -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  family <- .memwas_normalize_family(settings$family)
  approximation <- .memwas_normalize_approximation(settings$approximation %||% "laplace")
  init_approximation <- .memwas_normalize_init_approximation(settings$init_approximation %||% "variational_inference")
  se_method <- .memwas_normalize_se_method(settings$se_method %||% "hessian")
  method <- toupper(settings$method %||% "ML")
  if (!method %in% c("ML", "REML")) stop("MEMWAS: `method` must be either 'ML' or 'REML'.", call. = FALSE)
  random_cov <- match.arg(settings$random_cov %||% "diagonal", c("diagonal", "unstructured"))
  .memwas_check_scalar_nonnegative(settings$L1_penalty, "L1_penalty")
  .memwas_check_scalar_nonnegative(settings$L2_penalty, "L2_penalty")
  if (family != "gaussian" && method == "REML") method <- "ML"

  settings$family <- family
  settings$approximation <- approximation
  settings$init_approximation <- init_approximation
  settings$se_method <- se_method
  settings$method <- method
  settings$random_cov <- random_cov
  settings$dot_spec <- .memwas_validate_dot_settings(
    dot_predictors = settings$dot_predictors %||% NULL,
    dot_alternative = settings$dot_alternative %||% NULL,
    dot_threshold = settings$dot_threshold %||% 0,
    dot_alpha = settings$dot_alpha %||% 0.05
  )

  spline_variables <- settings$spline_variables %||% character(0L)
  spline_info <- settings$spline_info %||% list()
  augmented <- .memwas_add_stored_splines(settings$data, settings$formula, spline_variables, spline_info, engine = backend)

  if (verbose && length(spline_variables) > 0L) {
    message("MEMWAS: applying stored spline terms for: ", paste(spline_variables, collapse = ", "))
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 4: Fit the final mixed-effects model with the requested covariance -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  message("03. Fitting the final mixed-effects model with the requested covariance...")

  fit <- .memwas_fit_model(
    formula = augmented$formula,
    data = augmented$data,
    id = settings$id, time = settings$time,
    random = settings$random, family = family,
    autocor = settings$autocor,
    approximation = approximation,
    init_approximation = init_approximation,
    se_method = se_method,
    L1_penalty = settings$L1_penalty,
    L2_penalty = settings$L2_penalty,
    control = settings$control,
    method = method,
    random_cov = random_cov,
    engine = backend)

  fit$dot_spec <- settings$dot_spec
  fit$coefficient_table <- .memwas_apply_dot_to_coefficient_table(
    fit$coefficient_table, settings$dot_spec, term_map = fit$coefficient_term_map
  )
  fit$dot_unmatched_predictors <- .memwas_warn_unmatched_dot_predictors(
    settings$dot_spec, fit$coefficient_table$term, term_map = fit$coefficient_term_map
  )

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 5: Assemble and return the MEMWAS fit object -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  fit$call <- user_call
  fit$settings <- settings
  fit$formal_formula <- augmented$formula
  fit$spline_terms <- augmented$added_terms
  fit$spline_variables <- spline_variables
  fit$manual_overrides <- manual_overrides
  fit$methodological_note <- if (family == "gaussian") {
    if (backend == "cpp") {
      "Gaussian mixed model fitted by direct marginal likelihood/REML using registered C++ numerical helpers where available."
    } else {
      "Gaussian mixed model fitted by direct marginal likelihood/REML in base R."
    }
  } else {
    paste0("Non-Gaussian ", family, " model fitted with variational-inference initialization followed by the native-C++ ",
           .memwas_approximation_label(approximation),
           " final marginal approximation. Approximation diagnostics are stored in `approximation_diagnostics`.")
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 6: Screen GLMM assumptions for the final fitted model -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  fit$assumption_check_spec <- .memwas_assumption_spec_from_settings(settings, default = "None")
  fit$assumption_check_methods <- settings$assumption_check_methods %||% .memwas_assumption_spec_table(fit$assumption_check_spec)
  if (.memwas_assumption_is_empty_spec(fit$assumption_check_spec)) {
    fit$assumption_checks <- list(
      stage = "fitted_model",
      alpha = .memwas_assumption_alpha(settings),
      methods = fit$assumption_check_spec,
      summary = .memwas_assumption_empty_table(),
      tables = list(),
      note = "No GLMM assumption screening methods were requested."
    )
  } else {
    if (verbose) message("04. Screening GLMM assumptions for the final fitted model...")
    assumption_checks_try <- try(
      .memwas_run_assumption_checks(
        fit,
        settings = settings,
        spec = fit$assumption_check_spec,
        stage = "fitted_model"
      ),
      silent = TRUE
    )
    fit$assumption_checks <- if (inherits(assumption_checks_try, "try-error")) {
      .memwas_assumption_checks_unavailable(
        fit$assumption_check_spec,
        settings = settings,
        reason = paste("GLMM assumption screening failed:", as.character(assumption_checks_try)),
        stage = "fitted_model"
      )
    } else {
      assumption_checks_try
    }
  }

  class(fit) <- c("MEMWAS_fit", class(fit))
  fit
}

# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# S3 METHOD: print summary for MEMWAS_fit objects -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
