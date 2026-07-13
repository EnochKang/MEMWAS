#' @title Configure a MEMWAS mixed-effects analysis
#' @description Validate model inputs, optionally screen numeric fixed-effect predictors for non-linearity with restricted cubic splines, and return a `MEMWAS` settings object for `fit_MEMWAS()` or `tune_MEMWAS()`.
#' @param formula formula. Fixed-effect model formula such as `y ~ x1 + x2`; simple numeric fixed-effect terms may be screened for spline non-linearity.
#' @param family character or family. Response family. Supported options are `"gaussian"`, `"binomial"`, `"poisson"`, `"negative_binomial"`, `"gamma"`, and `"exponential"`; aliases include `"normal"`, `"logit"`, `"logistic"`, `"count"`, `"negbin"`, `"nb"`, and `"exp"`.
#' @param data data.frame. Long-format analysis data containing the response, predictors, subject identifier, time variable, and any random-effect variables.
#' @param id character. Name of the subject identifier column; grouped operations use this column to keep repeated measurements together.
#' @param time character. Name of the within-subject measurement-time column used for sorting and autocorrelation structures.
#' @param random formula or character. Random-effect design specification; `~ 1` fits a random intercept and strings such as `"~ 1 + time"` are accepted. A mixed-model style term like `~ 1 + time | id` is reduced to its left-hand random-effect design.
#' @param autocor character. One legacy within-subject residual correlation structure. Supported options include `"NONE"`, `"AR(1)"`, `"AR(p)"`, `"ARMA(1,1)"`, `"CS"`, `"TOEP"`, and `"UN"`. A fully named character vector may also be used as multiple-component shorthand.
#' @param serial Optional serial-process component or named list/vector of components created with `serial_component()`. Each component has its own one-column design and autocorrelation structure. Do not combine this with a non-`"NONE"` scalar `autocor` value.
#' @param L1_penalty numeric. Non-negative L1 fixed-effect penalty. The intercept is not penalized.
#' @param L2_penalty numeric. Non-negative L2 fixed-effect penalty. The intercept is not penalized.
#' @param control list. Advanced controls passed to fitting and screening, including optimizer controls, approximation controls, `continuous_time`, `toep_lags`, `max_unstructured_dim`, `spline_grid`, and `min_unique_nonlinear`. Family-specific entries include `negative_binomial_theta`, `gamma_shape`, `estimate_dispersion`, `aghq_nodes`, `aghq_max_dim`, `bootstrap_B`, and profile-likelihood controls.
#' @param method character. Likelihood method for Gaussian models: `"ML"` or `"REML"`. Non-Gaussian models are fitted with approximate ML-style working likelihoods and are forced to `"ML"`.
#' @param random_cov character. Random-effect covariance parameterization: `"diagonal"` or `"unstructured"`.
#' @param approximation character. Final non-Gaussian approximation strategy. The default is `"laplace"` after variational-inference initialization. Supported options are `"laplace"`, `"variational_inference"`, `"adaptive_gauss_hermite_quadrature"`, `"adaptive_gaussian_quadrature"`, `"saddlepoint"`, `"skew_corrected_laplace"`, and `"pql"`.
#' @param init_approximation character. Initial approximation for non-Gaussian models. Defaults to `"variational_inference"`; this is used for starting values before final Laplace/AGHQ/saddlepoint/skew-corrected fitting.
#' @param se_method character. Standard-error method for non-Gaussian final fits: `"hessian"`, `"cluster_sandwich"`, `"parametric_bootstrap"`, `"profile"`, or a character vector.
#' @param dot_predictors character. Optional fixed-effect coefficient terms for directional one-tailed tests (DOT). Names should match rows in the fixed-effect coefficient table, such as `"x1"` or `"groupB"`.
#' @param dot_alternative character, named character, list, or data.frame. DOT alternative(s): `"greater"` tests `beta > abs(dot_threshold)`, `"less"` tests `beta < -abs(dot_threshold)`, and `"two.sided"` uses the usual two-sided Wald test when `dot_threshold = 0` or a two-sided minimum-effect test of `|beta| > abs(dot_threshold)` when `dot_threshold > 0`. Named values override specific terms; a scalar value applies to `dot_predictors`, or to all non-intercept fixed-effect coefficient terms when `dot_predictors` is `NULL`. Data frames may be supplied to `dot_alternative`, `dot_threshold`, or `dot_alpha` with a `term`/`predictor` column plus relevant `alternative`, `threshold`/`dot_threshold`, or `alpha` columns.
#' @param dot_threshold numeric or data.frame. Predictor-specific minimally important value(s). Values are interpreted as magnitudes and converted with `abs()`: `"greater"` tests against `+abs(dot_threshold)`, `"less"` tests against `-abs(dot_threshold)`, and `"two.sided"` with a nonzero value tests whether the absolute effect exceeds that magnitude. Use a scalar, a vector aligned with `dot_predictors`, a named vector/list keyed by coefficient or model term, or a data frame with `term`/`predictor` and `threshold`/`dot_threshold` columns. Defaults to 0.
#' @param dot_alpha numeric. Significance level(s) for DOT. Use a scalar, a vector aligned with `dot_predictors`, or a named vector/list keyed by coefficient term. Defaults to 0.05.
#' @param nonlinear_alpha numeric. Significance threshold for likelihood-ratio screening of spline terms.
#' @param spline_knots integer. Number of restricted cubic spline knots; must be at least 3.
#' @param spline_probs numeric. Optional quantile probabilities used to place spline knots; length must equal `spline_knots`.
#' @param nonlinear_predictors character. Optional subset of fixed-effect predictor names to screen for non-linearity.
#' @param screen_nonlinear logical. Whether to perform spline screening during setup.
#' @param autocorrelation_check character or logical. Base-R residual autocorrelation screening method(s). Use \code{"All"}, \code{"None"}, \code{TRUE}, \code{FALSE}, or any combination of \code{"DurbinWatson"}, \code{"LjungBox"}, \code{"Lag1Correlation"}, and \code{"Runs"}. Default \code{"All"} runs all listed tests.
#' @param distribution_link_check character or logical. Base-R screening method(s) for the GLMM response distribution and link function. Use \code{"All"}, \code{"None"}, or any combination of \code{"PearsonDispersion"}, \code{"DevianceDispersion"}, \code{"LinkTest"}, \code{"QuantileResidualNormality"}, and \code{"GroupedCalibration"}. Default \code{"All"}.
#' @param conditional_independence_check character or logical. Base-R conditional residual-independence screening method(s). Use \code{"All"}, \code{"None"}, or any combination of \code{"WithinClusterLag1"}, \code{"ClusterMeanResidual"}, and \code{"Runs"}. Default \code{"All"}.
#' @param random_effects_normality_check character or logical. Base-R screening method(s) for empirical random-effect normality. Use \code{"All"}, \code{"None"}, or any combination of \code{"Shapiro"}, \code{"JarqueBera"}, and \code{"SkewKurtosis"}. Default \code{"All"}.
#' @param random_effects_predictor_independence_check character or logical. Base-R screening method(s) for independence between empirical random effects and predictors. Use \code{"All"}, \code{"None"}, or any combination of \code{"GroupMeanAssociation"}, \code{"Correlation"}, and \code{"RankCorrelation"}. Default \code{"All"}.
#' @param homogeneity_variance_check character or logical. Base-R screening method(s) for homogeneous variance / scale. Use \code{"All"}, \code{"None"}, or any combination of \code{"BreuschPagan"}, \code{"White"}, \code{"LeveneFitted"}, and \code{"LeveneGroup"}. Default \code{"All"}.
#' @param engine character. Execution engine: `"R"` runs all calculations in base R; `"cpp"` routes available numerical helpers, linear algebra, spline bases, penalties, metrics, and optimization kernels through registered C++ `.Call` routines with no Rcpp dependency.
#' @param verbose logical. Whether to print screening progress messages.
#' @param ... list. Additional user metadata stored in the returned settings object.
#' @details
#' `set_MEMWAS()` is the configuration layer for MEMWAS. It normalizes the response family, likelihood method, random-effect covariance option, and residual autocorrelation structure; checks that formula variables exist; optionally performs one-predictor-at-a-time restricted-cubic-spline screening; stores selected spline bases for the final formal fit, and performs dependency-free screening of six GLMM assumptions from the baseline mixed-effects screening model.
#'
#' Algorithm choices are controlled primarily by `family`, `method`, `random_cov`, `autocor`/`serial`, `approximation`, and `engine`. Gaussian models use marginal ML or REML. Binomial, Poisson, negative-binomial, gamma, and exponential models use variational-inference initialization followed by the selected native-C++ marginal approximation inside `fit_MEMWAS()`, so non-linearity screening for those families is approximate and is labeled by the selected final approximation method. Autocorrelation choices include independence, AR-type, compound-symmetry, Toeplitz, and unstructured residual correlation structures. When `engine = "cpp"`, numerical helper routines are dispatched to native C++ implementations registered with `.Call`; model orchestration, formula handling, and S3 object construction remain in R.
#' @returns MEMWAS. A settings object containing validated inputs, spline-screening results, selected spline metadata, GLMM assumption-screening specifications and results, and fitting controls.
#' @examples
#' \dontrun{
#' dat <- simulate_panel_data(n_id = 100, n_time = 5)
#' setup <- set_MEMWAS(var_1 ~ var_2 + var_3, data = dat, id = "id", time = "time",
#'                     family = "gaussian", autocor = "AR(1)", engine = "cpp")
#' fit <- fit_MEMWAS(setup)
#' }
#' @export
set_MEMWAS <- function(formula, family = "gaussian", data, id, time, random = ~ 1,
                       autocor = "AR(1)", L1_penalty = 0, L2_penalty = 0,
                       control = list(), method = "ML",
                       random_cov = c("diagonal", "unstructured"),
                       approximation = "laplace",
                       init_approximation = "variational_inference",
                       se_method = "hessian",
                       dot_predictors = NULL, dot_alternative = NULL,
                       dot_threshold = 0, dot_alpha = 0.05,
                       nonlinear_alpha = 0.05, spline_knots = 4L,
                       spline_probs = NULL, nonlinear_predictors = NULL,
                       screen_nonlinear = TRUE,
                       autocorrelation_check = "All",
                       distribution_link_check = "All",
                       conditional_independence_check = "All",
                       random_effects_normality_check = "All",
                       random_effects_predictor_independence_check = "All",
                       homogeneity_variance_check = "All",
                       engine = c("R", "cpp"),
                       serial = NULL,
                       verbose = TRUE, ...) {
  backend <- .memwas_normalize_engine(engine)

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 1: Validate user-supplied settings -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  message("01. Checking arguments...")

  if (missing(formula)) stop("Required argument `formula` is missing.", call. = FALSE)
  if (missing(data)) stop("Required argument `data` is missing.", call. = FALSE)
  if (missing(id)) stop("Required argument `id` is missing.", call. = FALSE)
  if (missing(time)) stop("Required argument `time` is missing.", call. = FALSE)
  if (!exists("fit_MEMWAS", mode = "function")) stop("Function `fit_MEMWAS()` must be defined before calling `set_MEMWAS()`.", call. = FALSE)
  if (!inherits(formula, "formula")) stop("`formula` must be an R formula.", call. = FALSE)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(id) || length(id) != 1L || !id %in% names(data)) {
    stop("`id` must be the name of an existing subject-id column in `data`.", call. = FALSE)
  }
  if (!is.character(time) || length(time) != 1L || !time %in% names(data)) {
    stop("`time` must be the name of an existing measurement-time column in `data`.", call. = FALSE)
  }
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  if (!is.null(serial)) {
    if (!exists(".memwas_validate_serial_syntax", mode = "function")) {
      stop("The multiple-serial engine is unavailable in this MEMWAS installation.", call. = FALSE)
    }
    .memwas_validate_serial_syntax(serial)
    if (missing(autocor)) {
      autocor <- "NONE"
    } else if (!.memwas_is_none_serial(autocor)) {
      stop("Supply all serial processes through `serial`; do not also supply a non-NONE scalar `autocor`.", call. = FALSE)
    }
  }
  if (!is.numeric(L1_penalty) || length(L1_penalty) != 1L || is.na(L1_penalty) || L1_penalty < 0) {
    stop("`L1_penalty` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(L2_penalty) || length(L2_penalty) != 1L || is.na(L2_penalty) || L2_penalty < 0) {
    stop("`L2_penalty` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(nonlinear_alpha) || length(nonlinear_alpha) != 1L || is.na(nonlinear_alpha) ||
      nonlinear_alpha <= 0 || nonlinear_alpha >= 1) {
    stop("`nonlinear_alpha` must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(spline_knots) || length(spline_knots) != 1L || is.na(spline_knots) || spline_knots < 3L) {
    stop("`spline_knots` must be at least 3.", call. = FALSE)
  }
  spline_knots <- as.integer(spline_knots)
  random_cov <- match.arg(random_cov)
  approximation <- .memwas_normalize_approximation(approximation)
  dot_spec <- .memwas_validate_dot_settings(dot_predictors = dot_predictors,
                                            dot_alternative = dot_alternative,
                                            dot_threshold = dot_threshold,
                                            dot_alpha = dot_alpha)
  assumption_check_spec <- .memwas_validate_assumption_check_settings(
    autocorrelation_check = autocorrelation_check,
    distribution_link_check = distribution_link_check,
    conditional_independence_check = conditional_independence_check,
    random_effects_normality_check = random_effects_normality_check,
    random_effects_predictor_independence_check = random_effects_predictor_independence_check,
    homogeneity_variance_check = homogeneity_variance_check
  )
  assumption_check_methods <- .memwas_assumption_spec_table(assumption_check_spec)

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 2: Use internal helper functions for setup, spline screening, and comparison -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  family <- .memwas_normalize_family(family)
  approximation <- .memwas_normalize_approximation(approximation)
  init_approximation <- .memwas_normalize_init_approximation(init_approximation)
  se_method <- .memwas_normalize_se_method(se_method)
  method <- toupper(method)
  if (!method %in% c("ML", "REML")) stop("`method` must be either 'ML' or 'REML'.", call. = FALSE)
  if (family != "gaussian" && method == "REML") {
    warning("REML is only used for Gaussian models; switching non-Gaussian formal analysis to ML with the selected approximation.", call. = FALSE)
    method <- "ML"
  }

  random_formula <- .memwas_parse_random_formula(random)
  missing_random_vars <- setdiff(all.vars(random_formula), names(data))
  if (length(missing_random_vars) > 0L) {
    stop("Random-effects variable(s) absent from `data`: ", paste(missing_random_vars, collapse = ", "),
         call. = FALSE)
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 3: Identify numeric fixed-effect predictors eligible for non-linear screening -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  message("02. Setting MEMWAS...")

  candidates <- .memwas_candidate_predictors(formula, data, id, time, nonlinear_predictors)
  if (verbose && screen_nonlinear) {
    if (length(candidates) == 0L) {
      message("set_MEMWAS: no simple numeric fixed-effect predictors were eligible for spline screening.")
    } else {
      message("set_MEMWAS: spline screening candidates: ", paste(candidates, collapse = ", "))
    }
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 4: Fit the baseline linear mixed-effects model for screening and Compare baseline against one-predictor-at-a-time spline models -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  nonlinear_summary <- data.frame()
  turning_points <- list()
  spline_info_all <- list()
  significant_splines <- character(0L)
  baseline_metrics <- NULL
  baseline_fit <- NULL

  if (isTRUE(screen_nonlinear) && length(candidates) > 0L) {

    message(" - Fitting the baseline linear mixed-effects model...")

    baseline_fit <- try(.memwas_fit_screening_model(formula, data, family = family, id = id, time = time, random = random, autocor = autocor, serial = serial, L1_penalty = L1_penalty, L2_penalty = L2_penalty, control = control, screen_method = "ML", random_cov = random_cov, approximation = approximation, init_approximation = init_approximation, se_method = se_method, engine = backend), silent = TRUE)
    if (inherits(baseline_fit, "try-error")) {
      stop("Baseline model failed during non-linear screening: ", as.character(baseline_fit), call. = FALSE)
    }
    baseline_metrics <- baseline_fit$metrics

    rows <- vector("list", length(candidates))
    names(rows) <- candidates
    min_unique <- control$min_unique_nonlinear %||% max(5L, spline_knots)
    grid_n <- max(25L, as.integer(control$spline_grid %||% 250L))

    for (v in candidates) {
      message(" - Iterating for checking non-linearity: ", v, "...")
      x <- data[[v]]
      if (length(unique(x[is.finite(x)])) < min_unique) {
        rows[[v]] <- .memwas_spline_failure_row(v, paste0("Skipped: fewer than ", min_unique, " unique finite values."), n = nrow(data))
        turning_points[[v]] <- numeric(0L)
        next
      }
      knots <- .memwas_get_spline_knots(x, spline_knots, spline_probs = spline_probs)
      if (length(knots) < 3L) {
        rows[[v]] <- .memwas_spline_failure_row(v, "Skipped: fewer than three unique spline knots after quantile calculation.", n = nrow(data))
        turning_points[[v]] <- numeric(0L)
        next
      }
      B <- .memwas_rcs_basis(x, knots, engine = backend)
      keep <- apply(B, 2L, function(z) stats::var(z, na.rm = TRUE) > 1e-12)
      B <- B[, keep, drop = FALSE]
      if (ncol(B) == 0L) {
        rows[[v]] <- .memwas_spline_failure_row(v, "Skipped: spline basis columns had near-zero variance.", n = nrow(data))
        turning_points[[v]] <- numeric(0L)
        next
      }
      basis_names <- .memwas_make_basis_names(v, ncol(B), names(data))
      data_s <- data
      for (j in seq_len(ncol(B))) data_s[[basis_names[j]]] <- B[, j]
      formula_s <- .memwas_add_terms_to_formula(formula, basis_names)
      fit_s <- try(.memwas_fit_screening_model(formula_s, data_s, family = family, id = id, time = time, random = random, autocor = autocor, serial = serial, L1_penalty = L1_penalty, L2_penalty = L2_penalty, control = control, screen_method = "ML", random_cov = random_cov, approximation = approximation, init_approximation = init_approximation, se_method = se_method, engine = backend), silent = TRUE)
      if (inherits(fit_s, "try-error")) {
        rows[[v]] <- .memwas_spline_failure_row(v, paste0("Spline model failed: ", as.character(fit_s)), n = nrow(data))
        turning_points[[v]] <- numeric(0L)
        next
      }
      ll0 <- baseline_fit$metrics$logLik
      ll1 <- fit_s$metrics$logLik
      dflrt <- ncol(B)
      lrt <- 2 * (ll1 - ll0)
      pval <- if (is.finite(lrt)) stats::pchisq(max(lrt, 0), df = dflrt, lower.tail = FALSE) else NA_real_
      is_nl <- is.finite(pval) && pval < nonlinear_alpha
      tp <- .memwas_spline_turning_points(v, fit_s, knots, basis_names, data, grid_n, engine = backend)
      turning_points[[v]] <- tp
      spline_info_all[[v]] <- list(knots = knots, basis_names = basis_names,
                                   df = ncol(B), screening_p_value = pval,
                                   turning_points = tp)
      if (is_nl) significant_splines <- c(significant_splines, v)
      rows[[v]] <- data.frame(
        predictor = v,
        n = nrow(data),
        knots = paste(signif(knots, 6), collapse = "; "),
        linear_logLik = as.numeric(ll0),
        spline_logLik = as.numeric(ll1),
        linear_AIC = as.numeric(baseline_fit$metrics$AIC),
        spline_AIC = as.numeric(fit_s$metrics$AIC),
        linear_BIC = as.numeric(baseline_fit$metrics$BIC),
        spline_BIC = as.numeric(fit_s$metrics$BIC),
        delta_logLik = as.numeric(ll1 - ll0),
        delta_AIC = as.numeric(fit_s$metrics$AIC - baseline_fit$metrics$AIC),
        delta_BIC = as.numeric(fit_s$metrics$BIC - baseline_fit$metrics$BIC),
        LRT = as.numeric(lrt),
        df = as.integer(dflrt),
        p_value = as.numeric(pval),
        nonlinear = isTRUE(is_nl),
        turning_points = if (length(tp)) paste(signif(tp, 6), collapse = "; ") else "",
        convergence_linear = as.integer(baseline_fit$convergence),
        convergence_spline = as.integer(fit_s$convergence),
        message = if (fit_s$approximate) paste0("Approximate ", fit_s$approximation_label %||% fit_s$approximation %||% "non-Gaussian", " comparison for non-Gaussian model.") else "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    nonlinear_summary <- do.call(rbind, rows)
    row.names(nonlinear_summary) <- NULL
  } else {
    baseline_fit <- try(.memwas_fit_screening_model(formula, data, family = family, id = id, time = time, random = random, autocor = autocor, serial = serial, L1_penalty = L1_penalty, L2_penalty = L2_penalty, control = control, screen_method = "ML", random_cov = random_cov, approximation = approximation, init_approximation = init_approximation, se_method = se_method, engine = backend), silent = TRUE)
    if (!inherits(baseline_fit, "try-error")) baseline_metrics <- baseline_fit$metrics
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 5: Screen GLMM assumptions on the baseline mixed-effects model -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  assumption_checks <- .memwas_assumption_checks_unavailable(
    assumption_check_spec,
    settings = list(control = control),
    reason = "Baseline mixed-effects model was unavailable, so GLMM assumption screening was not run.",
    stage = "setup_baseline_model"
  )
  if (.memwas_assumption_is_empty_spec(assumption_check_spec)) {
    assumption_checks <- list(
      stage = "setup_baseline_model",
      alpha = .memwas_assumption_alpha(list(control = control)),
      methods = assumption_check_spec,
      summary = .memwas_assumption_empty_table(),
      tables = list(),
      note = "No GLMM assumption screening methods were requested."
    )
  } else if (!is.null(baseline_fit) && !inherits(baseline_fit, "try-error") && inherits(baseline_fit, "MEMWAS_fit")) {
    if (verbose) message(" - Screening GLMM assumptions on the baseline mixed-effects model...")
    settings_for_checks <- baseline_fit$settings %||% list(control = control, data = data, id = id, time = time, family = family)
    assumption_checks_try <- try(
      .memwas_run_assumption_checks(
        baseline_fit,
        settings = settings_for_checks,
        spec = assumption_check_spec,
        stage = "setup_baseline_model"
      ),
      silent = TRUE
    )
    if (inherits(assumption_checks_try, "try-error")) {
      assumption_checks <- .memwas_assumption_checks_unavailable(
        assumption_check_spec,
        settings = settings_for_checks,
        reason = paste("GLMM assumption screening failed:", as.character(assumption_checks_try)),
        stage = "setup_baseline_model"
      )
    } else {
      assumption_checks <- assumption_checks_try
    }
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 6: Store significant non-linear predictors for the formal MEMWAS fit -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  significant_splines <- unique(significant_splines)
  spline_info <- spline_info_all[significant_splines]
  formal_formula <- formula
  if (length(significant_splines) > 0L) {
    formal_terms <- unlist(lapply(significant_splines, function(v) spline_info[[v]]$basis_names), use.names = FALSE)
    formal_formula <- .memwas_add_terms_to_formula(formula, formal_terms)
    if (verbose) message("set_MEMWAS: formal MEMWAS will include spline terms for: ", paste(significant_splines, collapse = ", "))
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # Step 7: Return the MEMWAS class object -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

  out <- list(
    call = match.call(),
    formula = formula,
    formal_formula = formal_formula,
    family = family,
    data = data,
    id = id,
    time = time,
    random = random,
    random_cov = random_cov,
    approximation = approximation,
    init_approximation = init_approximation,
    se_method = se_method,
    engine = backend,
    autocor = autocor,
    serial = serial,
    L1_penalty = L1_penalty,
    L2_penalty = L2_penalty,
    control = control,
    method = method,
    dot_predictors = dot_predictors,
    dot_alternative = dot_alternative,
    dot_threshold = dot_threshold,
    dot_alpha = dot_alpha,
    dot_spec = dot_spec,
    nonlinear_alpha = nonlinear_alpha,
    spline_knots = spline_knots,
    spline_probs = spline_probs,
    nonlinear_predictors = nonlinear_predictors,
    screen_nonlinear = screen_nonlinear,
    autocorrelation_check = autocorrelation_check,
    distribution_link_check = distribution_link_check,
    conditional_independence_check = conditional_independence_check,
    random_effects_normality_check = random_effects_normality_check,
    random_effects_predictor_independence_check = random_effects_predictor_independence_check,
    homogeneity_variance_check = homogeneity_variance_check,
    assumption_check_spec = assumption_check_spec,
    assumption_check_methods = assumption_check_methods,
    assumption_checks = assumption_checks,
    assumption_check_model = "setup_baseline_model",
    nonlinear_summary = nonlinear_summary,
    turning_points = turning_points,
    spline_variables = significant_splines,
    spline_info = spline_info,
    all_screened_spline_info = spline_info_all,
    baseline_screen_metrics = baseline_metrics,
    extra = list(...),
    verbose = verbose,
    note = paste(
      "Gaussian models use direct marginal ML/REML.",
      "Non-Gaussian models use variational-inference initialization followed by the selected native-C++ marginal approximation; Laplace is the default final approximation."
    )
  )
  class(out) <- "MEMWAS"
  out
}


# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# S3 METHOD: print summary for set_MEMWAS objects  -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
