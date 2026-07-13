#' @title Tune MEMWAS elastic-net penalties by grouped K-fold cross-validation
#' @description Tune MEMWAS fixed-effect L1/L2 penalties using subject-grouped K-fold cross-validation and a lightweight surrogate search over `lambda` and `alpha`. The returned tuning object now stores a structured `best_model` summary with final fixed-effect coefficients, random-effect BLUPs and covariance estimates, penalty formulas and values, and autocorrelation formulas and penalty values when the final model is refitted.
#' @param object MEMWAS or formula or NULL. A settings object from `set_MEMWAS()`, a formula, or `NULL` when direct model arguments are supplied.
#' @param formula formula. Fixed-effect model formula used when `object` is not a populated settings object.
#' @param family character or family. Response family: `"gaussian"`, `"binomial"`, `"poisson"`, `"negative_binomial"`, `"gamma"`, or `"exponential"`.
#' @param data data.frame. Long-format analysis data.
#' @param id character. Subject identifier column used for grouped cross-validation folds.
#' @param time character. Measurement-time column.
#' @param random formula or character. Random-effect design.
#' @param autocor character. Legacy residual autocorrelation structure to use in each cross-validation fit, or a fully named multiple-component shorthand.
#' @param serial Optional serial-process component or named list/vector of components created with `serial_component()`. The same definitions are retained in every grouped fold and the final refit.
#' @param control list. Fitting, approximation, optimizer, spline, and autocorrelation controls passed to `fit_MEMWAS()`.
#' @param method character. `"ML"` or `"REML"` for Gaussian models; non-Gaussian models are evaluated with ML and the selected approximation.
#' @param random_cov character. Random-effect covariance option: `"diagonal"` or `"unstructured"`.
#' @param approximation character. Final non-Gaussian approximation strategy passed to each fold fit and final refit. The default is `"laplace"` after variational-inference initialization. Supported options are `"laplace"`, `"adaptive_gauss_hermite_quadrature"`, `"adaptive_gaussian_quadrature"`, `"saddlepoint"`, `"skew_corrected_laplace"`, `"variational_inference"`, and legacy `"pql"`.
#' @param init_approximation character. Initial approximation used for non-Gaussian fold fits and final refit. Defaults to `"variational_inference"`.
#' @param se_method character. Standard-error method for non-Gaussian final refits. Use `"hessian"`, `"cluster_sandwich"`, `"parametric_bootstrap"`, `"profile"`, or a character vector.
#' @param dot_predictors character. Optional fixed-effect coefficient terms for directional one-tailed tests (DOT) in the final refit and any stored fold fits. Names should match rows in the fixed-effect coefficient table.
#' @param dot_alternative character, named character, list, or data.frame. DOT alternative(s): `"greater"` tests `beta > abs(dot_threshold)`, `"less"` tests `beta < -abs(dot_threshold)`, and `"two.sided"` uses the usual two-sided Wald test when `dot_threshold = 0` or a two-sided minimum-effect test of `|beta| > abs(dot_threshold)` when `dot_threshold > 0`. Named values override specific terms; a scalar value applies to `dot_predictors`, or to all non-intercept fixed-effect coefficient terms when `dot_predictors` is `NULL`. Data frames may be supplied to `dot_alternative`, `dot_threshold`, or `dot_alpha` with a `term`/`predictor` column plus relevant `alternative`, `threshold`/`dot_threshold`, or `alpha` columns.
#' @param dot_threshold numeric or data.frame. Predictor-specific minimally important value(s). Values are interpreted as magnitudes and converted with `abs()`: `"greater"` tests against `+abs(dot_threshold)`, `"less"` tests against `-abs(dot_threshold)`, and `"two.sided"` with a nonzero value tests whether the absolute effect exceeds that magnitude. Use a scalar, a vector aligned with `dot_predictors`, a named vector/list keyed by coefficient or model term, or a data frame with `term`/`predictor` and `threshold`/`dot_threshold` columns. Defaults to 0.
#' @param dot_alpha numeric. Significance level(s) for DOT. Use a scalar, a vector aligned with `dot_predictors`, or a named vector/list keyed by coefficient term. Defaults to 0.05.
#' @param engine character. Execution engine: `"R"` for base-R numerical helpers or `"cpp"` for registered native C++ helper routines in every cross-validation fit and tuning kernel where available.
#' @param K integer. Number of grouped cross-validation folds; all rows for a subject stay in one fold.
#' @param metric character. Prediction error metric to minimize: `"MAE"`, `"MSE"`, `"RMSE"`, `"MAPE"`, or `"SMAPE"`.
#' @param stability_metric character. Dispersion metric for fold stability: `"IQR"`, `"SD"`, or `"variance"`.
#' @param surrogate character. Surrogate model for proposing new penalty values: `"QRS"` for quadratic response surface or `"GP"` for a small Gaussian-process surrogate.
#' @param use_stability_penalty logical. Whether to add a stability penalty to the median cross-validation error.
#' @param lambda_log_range numeric. Two-element range for `log(lambda, base = lambda_base)`.
#' @param lambda_base numeric. Base used to transform log-lambda values to lambda.
#' @param alpha_range numeric. Two-element range for elastic-net mixing parameter `alpha`, constrained to `[0, 1]`.
#' @param initial_design character. Initial design strategy: `"lhs"` for Latin-hypercube-like sampling plus boundaries or `"boundary"` for boundary points only.
#' @param n_initial integer. Number of initial penalty pairs to evaluate for the `"lhs"` design.
#' @param initial_points data.frame. Optional user-supplied starting points with `log_lambda` or `lambda`, and `alpha`.
#' @param max_iter integer. Maximum number of surrogate-guided search iterations after initial evaluations.
#' @param improvement_tol numeric. Minimum improvement in selection score considered meaningful.
#' @param patience integer. Stop after this many consecutive negligible-improvement iterations.
#' @param candidate_grid_size integer. Number of candidate penalty pairs sampled at each surrogate iteration.
#' @param exploration numeric. Exploration multiplier for surrogate uncertainty.
#' @param stability_weight numeric. Weight for the stability penalty in the selection score.
#' @param seed integer. Optional random seed for reproducible fold assignment and design generation.
#' @param refit_final logical. Whether to refit the best MEMWAS model on the full data set.
#' @param autocor_regularization numeric or list. Optional regularization applied to autocorrelation parameters during fitting.
#' @param autocor_regularization_type character. Autocorrelation regularization type: `"L2"`, `"L1"`, `"elasticnet"`, or `"none"`.
#' @param autocor_regularization_alpha numeric. Mixing value for elastic-net autocorrelation regularization.
#' @param mape_epsilon numeric. Positive stabilizer used in MAPE and SMAPE denominators.
#' @param fail_value numeric. Metric value assigned to failed folds.
#' @param keep_fold_fits logical. Whether to store all fold-level fitted model objects.
#' @param verbose logical. Whether to print tuning progress.
#' @param ... list. Additional arguments passed to `fit_MEMWAS()`.
#' @details
#' `tune_MEMWAS()` tunes fixed-effect penalties through the reparameterization `L1_penalty = alpha * lambda` and `L2_penalty = (1 - alpha) * lambda`. Cross-validation is grouped by subject ID, so no subject contributes rows to both training and validation within the same fold. The best row minimizes median validation error plus an optional stability penalty computed from fold-level error dispersion.
#'
#' Search begins with boundary and optional Latin-hypercube-like initial points. It then fits either a quadratic response-surface surrogate (`surrogate = "QRS"`) or a compact Gaussian-process surrogate (`surrogate = "GP"`) over observed cross-validation results and proposes new penalty pairs with uncertainty-guided exploration. `metric`, `stability_metric`, `surrogate`, `lambda_log_range`, `alpha_range`, and `engine` are the main tuning algorithm controls.
#'
#' With `engine = "cpp"`, available numerical work inside each MEMWAS fit and several tuning kernels, including metrics, fold assignment, and covariance-kernel calculations, are dispatched to compiled C++ routines through `.Call` without using Rcpp or other external packages.
#' @returns tune_MEMWAS. A tuning object containing the best penalty pair, all observed results, fold-level cross-validation details, optional final fit, a structured `best_model` component with final formulas and coefficients, fold assignments, settings, tuning controls, and convergence diagnostics.
#' @examples
#' \dontrun{
#' dat <- simulate_panel_data(n_id = 50, n_time = 5)
#' setup <- set_MEMWAS(var_1 ~ var_2 + var_3, data = dat, id = "id", time = "time",
#'                     screen_nonlinear = FALSE)
#' tuned <- tune_MEMWAS(setup, K = 3, max_iter = 2, n_initial = 5, engine = "cpp")
#' tuned$best
#' }
#' @export
tune_MEMWAS <- function(object = NULL,
                        formula = NULL,
                        family = NULL,
                        data = NULL,
                        id = NULL,
                        time = NULL,
                        random = NULL,
                        autocor = NULL,
                        control = NULL,
                        method = NULL,
                        random_cov = NULL,
                        approximation = NULL,
                        init_approximation = NULL,
                        se_method = NULL,
                        dot_predictors = NULL,
                        dot_alternative = NULL,
                        dot_threshold = 0,
                        dot_alpha = 0.05,
                        engine = NULL,
                        K = 5L,
                        metric = c("MAE", "MSE", "RMSE", "MAPE", "SMAPE"),
                        stability_metric = c("IQR", "SD", "variance"),
                        surrogate = c("QRS", "GP"),
                        use_stability_penalty = TRUE,
                        lambda_log_range = c(-6, 2),
                        lambda_base = 10,
                        alpha_range = c(0, 1),
                        initial_design = c("lhs", "boundary"),
                        n_initial = 10L,
                        initial_points = NULL,
                        max_iter = 20L,
                        improvement_tol = 1e-6,
                        patience = 3L,
                        candidate_grid_size = 1000L,
                        exploration = 1,
                        stability_weight = 1,
                        seed = NULL,
                        refit_final = TRUE,
                        autocor_regularization = NULL,
                        autocor_regularization_type = c("L2", "L1", "elasticnet", "none"),
                        autocor_regularization_alpha = 0,
                        mape_epsilon = .Machine$double.eps,
                        fail_value = Inf,
                        keep_fold_fits = FALSE,
                        serial = NULL,
                        verbose = TRUE,
                        ...) {
  call <- match.call()
  dots <- list(...)

  if (verbose) {
    message("01. Checking inputs and arguments...")
  }

  if (!exists("fit_MEMWAS", mode = "function")) {
    stop("Function `fit_MEMWAS()` must be defined before calling `tune_MEMWAS()`.", call. = FALSE)
  }

  metric <- toupper(as.character(metric)[1L])
  if (!metric %in% c("MAE", "MSE", "RMSE", "MAPE", "SMAPE")) {
    stop("`metric` must be one of 'MAE', 'MSE', 'RMSE', 'MAPE', or 'SMAPE'.", call. = FALSE)
  }
  stability_metric <- toupper(as.character(stability_metric)[1L])
  if (stability_metric == "VAR") stability_metric <- "VARIANCE"
  if (!stability_metric %in% c("IQR", "SD", "VARIANCE")) {
    stop("`stability_metric` must be one of 'IQR', 'SD', or 'variance'.", call. = FALSE)
  }
  surrogate <- toupper(as.character(surrogate)[1L])
  if (!surrogate %in% c("QRS", "GP")) {
    stop("`surrogate` must be either 'QRS' or 'GP'.", call. = FALSE)
  }
  initial_design <- tolower(as.character(initial_design)[1L])
  if (!initial_design %in% c("lhs", "boundary")) {
    stop("`initial_design` must be either 'lhs' or 'boundary'.", call. = FALSE)
  }
  autocor_regularization_type <- tolower(as.character(autocor_regularization_type)[1L])
  if (!autocor_regularization_type %in% c("l2", "l1", "elasticnet", "none")) {
    stop("`autocor_regularization_type` must be one of 'L2', 'L1', 'elasticnet', or 'none'.", call. = FALSE)
  }

  K <- as.integer(K[1L])
  n_initial <- as.integer(n_initial[1L])
  max_iter <- as.integer(max_iter[1L])
  patience <- as.integer(patience[1L])
  candidate_grid_size <- as.integer(candidate_grid_size[1L])

  if (!is.finite(K) || K < 2L) stop("`K` must be an integer of at least 2.", call. = FALSE)
  if (!is.finite(n_initial) || n_initial < 1L) stop("`n_initial` must be a positive integer.", call. = FALSE)
  if (!is.finite(max_iter) || max_iter < 0L) stop("`max_iter` must be a non-negative integer.", call. = FALSE)
  if (!is.finite(patience) || patience < 1L) stop("`patience` must be a positive integer.", call. = FALSE)
  if (!is.finite(candidate_grid_size) || candidate_grid_size < 10L) {
    stop("`candidate_grid_size` must be at least 10.", call. = FALSE)
  }
  if (!is.numeric(lambda_log_range) || length(lambda_log_range) != 2L || any(!is.finite(lambda_log_range))) {
    stop("`lambda_log_range` must contain two finite numeric values.", call. = FALSE)
  }
  if (!is.numeric(lambda_base) || length(lambda_base) != 1L || !is.finite(lambda_base) || lambda_base <= 1) {
    stop("`lambda_base` must be a single numeric value greater than 1.", call. = FALSE)
  }
  if (!is.numeric(alpha_range) || length(alpha_range) != 2L || any(!is.finite(alpha_range))) {
    stop("`alpha_range` must contain two finite numeric values.", call. = FALSE)
  }
  lambda_log_range <- sort(as.numeric(lambda_log_range))
  alpha_range <- sort(as.numeric(alpha_range))
  if (alpha_range[1L] < 0 || alpha_range[2L] > 1 || alpha_range[1L] == alpha_range[2L]) {
    stop("`alpha_range` must lie within [0, 1] and contain two distinct values.", call. = FALSE)
  }
  if (!is.numeric(exploration) || length(exploration) != 1L || !is.finite(exploration) || exploration < 0) {
    stop("`exploration` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(stability_weight) || length(stability_weight) != 1L || !is.finite(stability_weight) || stability_weight < 0) {
    stop("`stability_weight` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(improvement_tol) || length(improvement_tol) != 1L || !is.finite(improvement_tol) || improvement_tol < 0) {
    stop("`improvement_tol` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(mape_epsilon) || length(mape_epsilon) != 1L || !is.finite(mape_epsilon) || mape_epsilon <= 0) {
    stop("`mape_epsilon` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  supplied <- list(
    formula = !missing(formula), family = !missing(family), data = !missing(data),
    id = !missing(id), time = !missing(time), random = !missing(random),
    autocor = !missing(autocor), serial = !missing(serial),
    control = !missing(control), method = !missing(method),
    random_cov = !missing(random_cov), approximation = !missing(approximation), init_approximation = !missing(init_approximation), se_method = !missing(se_method), dot_predictors = !missing(dot_predictors),
    dot_alternative = !missing(dot_alternative), dot_threshold = !missing(dot_threshold),
    dot_alpha = !missing(dot_alpha), engine = !missing(engine)
  )

  if (missing(object) || is.null(object)) {
    settings <- .memwas_make_default_settings(call, verbose = FALSE)
  } else if (inherits(object, "MEMWAS")) {
    settings <- object
  } else if (inherits(object, "formula")) {
    settings <- .memwas_make_default_settings(call, verbose = FALSE)
    if (supplied$formula && !is.null(formula)) {
      stop("Supply the model formula either as `object` or `formula`, not both.", call. = FALSE)
    }
    settings$formula <- object
  } else {
    stop("`object` must be a MEMWAS settings object, a model formula, or NULL.", call. = FALSE)
  }

  if (supplied$formula) settings$formula <- formula
  if (supplied$family) settings$family <- family
  if (supplied$data) settings$data <- data
  if (supplied$id) settings$id <- id
  if (supplied$time) settings$time <- time
  if (supplied$random) settings$random <- random
  if (supplied$autocor) settings$autocor <- autocor
  if (supplied$serial) settings$serial <- serial
  if (!is.null(settings$serial %||% NULL)) {
    if (!exists(".memwas_validate_serial_syntax", mode = "function")) {
      stop("The multiple-serial engine is unavailable in this MEMWAS installation.", call. = FALSE)
    }
    .memwas_validate_serial_syntax(settings$serial)
    if (!supplied$autocor) {
      settings$autocor <- "NONE"
    } else if (!.memwas_is_none_serial(settings$autocor)) {
      stop("Supply all serial processes through `serial`; do not also supply a non-NONE scalar `autocor`.", call. = FALSE)
    }
  }
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
  if (supplied$engine) settings$engine <- engine

  if (is.null(settings$control)) settings$control <- list()
  if (is.null(settings$dot_threshold)) settings$dot_threshold <- 0
  if (is.null(settings$dot_alpha)) settings$dot_alpha <- 0.05
  if (!is.list(settings$control)) stop("`control` must be a list.", call. = FALSE)
  if (is.null(settings$random_cov)) settings$random_cov <- "diagonal"
  settings$random_cov <- match.arg(settings$random_cov, c("diagonal", "unstructured"))
  if (is.null(settings$method)) settings$method <- "ML"
  if (is.null(settings$engine)) settings$engine <- "R"
  settings$engine <- .memwas_normalize_engine(settings$engine)
  backend <- settings$engine
  settings$method <- toupper(settings$method)
  if (!settings$method %in% c("ML", "REML")) stop("`method` must be either 'ML' or 'REML'.", call. = FALSE)
  if (is.null(settings$family)) settings$family <- "gaussian"
  settings$family <- .memwas_normalize_family(settings$family)
  settings$approximation <- .memwas_normalize_approximation(settings$approximation %||% "laplace")
  settings$init_approximation <- .memwas_normalize_init_approximation(settings$init_approximation %||% "variational_inference")
  settings$se_method <- .memwas_normalize_se_method(settings$se_method %||% "hessian")
  if (settings$family != "gaussian" && settings$method == "REML") settings$method <- "ML"
  settings$dot_spec <- .memwas_validate_dot_settings(
    dot_predictors = settings$dot_predictors %||% NULL,
    dot_alternative = settings$dot_alternative %||% NULL,
    dot_threshold = settings$dot_threshold %||% 0,
    dot_alpha = settings$dot_alpha %||% 0.05
  )

  if (!inherits(settings$formula, "formula")) stop("A model `formula` must be supplied.", call. = FALSE)
  if (!is.data.frame(settings$data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(settings$id) || length(settings$id) != 1L || !settings$id %in% names(settings$data)) {
    stop("`id` must be the name of an existing subject-id column in `data`.", call. = FALSE)
  }
  if (!is.character(settings$time) || length(settings$time) != 1L || !settings$time %in% names(settings$data)) {
    stop("`time` must be the name of an existing measurement-time column in `data`.", call. = FALSE)
  }
  if (any(is.na(settings$data[[settings$id]]))) {
    stop("The subject-id column contains missing values; grouped K-fold splitting requires non-missing IDs.", call. = FALSE)
  }

  settings$control <- .memwas_augment_control(
    settings$control,
    autocor_regularization = autocor_regularization,
    autocor_regularization_type = autocor_regularization_type,
    autocor_regularization_alpha = autocor_regularization_alpha
  )

  fold_split <- .memwas_make_id_fold_split(settings$data[[settings$id]], K = K, engine = backend)
  row_folds <- fold_split$row_folds
  id_folds <- fold_split$id_folds
  fold_table <- table(row_folds)
  id_fold_table <- table(id_folds$fold)
  if (verbose) {
    message("tune_MEMWAS: Created ", K," grouped folds using unique `", settings$id, "` values; fold row counts: ",
            paste(as.integer(fold_table), collapse = ", "),
            "; fold subject counts: ",
            paste(as.integer(id_fold_table), collapse = ", "), ".\n")
  }

  initial <- .memwas_make_initial_points(lambda_log_range, alpha_range, initial_design = initial_design, n_initial = n_initial, initial_points = initial_points, lambda_base = lambda_base)
  observed <- data.frame()
  cv_detail_list <- list()
  fold_fit_list <- list()
  eval_counter <- 0L

  if (verbose) {
    message("02. Setting ", nrow(initial), " initial hyperparameter pair(s)...\n")
  }

  for (i in seq_len(nrow(initial))) {
    eval_counter <- eval_counter + 1L
    ev <- .memwas_evaluate_penalty_pair(initial$log_lambda[i], initial$alpha[i], iteration = 0L, stage = "initial",
                                      settings = settings, row_folds = row_folds, K = K,
                                      metric = metric, stability_metric = stability_metric,
                                      lambda_base = lambda_base, alpha_range = alpha_range,
                                      use_stability_penalty = use_stability_penalty,
                                      stability_weight = stability_weight, mape_epsilon = mape_epsilon,
                                      fail_value = fail_value, keep_fold_fits = keep_fold_fits,
                                      engine = backend, autocor_regularization = autocor_regularization,
                                      autocor_regularization_type = autocor_regularization_type,
                                      autocor_regularization_alpha = autocor_regularization_alpha, ...)
    ev$row$evaluation <- eval_counter
    ev$folds$evaluation <- eval_counter
    observed <- rbind(observed, ev$row)
    cv_detail_list[[length(cv_detail_list) + 1L]] <- ev$folds
    if (isTRUE(keep_fold_fits)) fold_fit_list[[as.character(eval_counter)]] <- ev$fold_fits
    if (verbose) {
      message("tune_MEMWAS: Initial ", i, " / ", nrow(initial), ": lambda = ", signif(ev$row$lambda, 4),
              ", alpha = ", signif(ev$row$alpha, 4),
              ", median ", metric, " = ", signif(ev$row$median_metric, 5),
              ", stability = ", signif(ev$row$stability, 5))
    }
  }

  finite_scores <- observed$selection_score[is.finite(observed$selection_score)]
  if (!length(finite_scores)) {
    stop("All initial MEMWAS cross-validation evaluations failed; inspect fold errors in the returned diagnostics or simplify the model.", call. = FALSE)
  }
  best_score <- min(finite_scores)
  no_improve <- 0L
  stop_reason <- "maximum iterations reached"
  completed_iter <- 0L

  for (iter in seq_len(max_iter)) {

    if (iter == 1 && verbose) {
      message(" ")
      message("03. Iterating for tuning...")
      message("Iteration ", iter, " / ", max_iter, " --------------------")
    } else {
      message(" ")
      message("Iteration ", iter, " / ", max_iter, " --------------------")
    }

    completed_iter <- iter
    candidates <- .memwas_not_observed(.memwas_make_candidate_points(candidate_grid_size, lambda_log_range, alpha_range), observed)
    if (nrow(candidates) == 0L) {
      stop_reason <- "no unevaluated candidate points remained"
      break
    }

    ok <- is.finite(observed$median_metric) & is.finite(observed$stability)
    Xobs <- as.matrix(observed[ok, c("log_lambda", "alpha"), drop = FALSE])
    yobs <- observed$median_metric[ok]
    sobs <- observed$stability[ok]

    if (nrow(Xobs) >= 2L) {
      perf_model <- .memwas_fit_surrogate(Xobs, yobs, surrogate = surrogate, lambda_log_range = lambda_log_range, alpha_range = alpha_range, engine = backend)
      Xcand <- as.matrix(candidates[, c("log_lambda", "alpha"), drop = FALSE])
      perf_pred <- .memwas_predict_surrogate(perf_model, Xcand, lambda_log_range = lambda_log_range, alpha_range = alpha_range, engine = backend)
      score <- perf_pred$mean - exploration * perf_pred$se
      stability_pred <- NULL
      if (isTRUE(use_stability_penalty)) {
        stab_model <- .memwas_fit_surrogate(Xobs, sobs, surrogate = surrogate, lambda_log_range = lambda_log_range, alpha_range = alpha_range, engine = backend)
        stability_pred <- .memwas_predict_surrogate(stab_model, Xcand, lambda_log_range = lambda_log_range, alpha_range = alpha_range, engine = backend)
        score <- score + stability_weight * pmax(stability_pred$mean, 0)
      }
      score[!is.finite(score)] <- Inf
      next_idx <- which.min(score)
      next_point <- candidates[next_idx, , drop = FALSE]
      next_stage <- "surrogate_ucb"
    } else {
      next_idx <- sample.int(nrow(candidates), 1L)
      next_point <- candidates[next_idx, , drop = FALSE]
      next_stage <- "random_fallback"
    }

    eval_counter <- eval_counter + 1L
    ev <- .memwas_evaluate_penalty_pair(next_point$log_lambda[1L], next_point$alpha[1L], iteration = iter, stage = next_stage,
                                      settings = settings, row_folds = row_folds, K = K,
                                      metric = metric, stability_metric = stability_metric,
                                      lambda_base = lambda_base, alpha_range = alpha_range,
                                      use_stability_penalty = use_stability_penalty,
                                      stability_weight = stability_weight, mape_epsilon = mape_epsilon,
                                      fail_value = fail_value, keep_fold_fits = keep_fold_fits,
                                      engine = backend, autocor_regularization = autocor_regularization,
                                      autocor_regularization_type = autocor_regularization_type,
                                      autocor_regularization_alpha = autocor_regularization_alpha, ...)
    ev$row$evaluation <- eval_counter
    ev$folds$evaluation <- eval_counter
    observed <- rbind(observed, ev$row)
    cv_detail_list[[length(cv_detail_list) + 1L]] <- ev$folds
    if (isTRUE(keep_fold_fits)) fold_fit_list[[as.character(eval_counter)]] <- ev$fold_fits

    new_finite_scores <- observed$selection_score[is.finite(observed$selection_score)]
    new_best <- if (length(new_finite_scores)) min(new_finite_scores) else best_score
    improvement <- best_score - new_best
    if (!is.finite(improvement) || improvement <= improvement_tol) {
      no_improve <- no_improve + 1L
    } else {
      no_improve <- 0L
      best_score <- new_best
    }

    if (verbose) {
      message("  iter ", iter, "/", max_iter, ": lambda=", signif(ev$row$lambda, 4),
              ", alpha=", signif(ev$row$alpha, 4),
              ", median ", metric, "=", signif(ev$row$median_metric, 5),
              ", stability=", signif(ev$row$stability, 5),
              ", selection score=", signif(ev$row$selection_score, 5))
    }

    if (no_improve >= patience) {
      stop_reason <- paste0("negligible improvement for ", patience, " consecutive iteration(s)")
      break
    }
  }

  observed <- observed[order(observed$selection_score, observed$median_metric, observed$stability), , drop = FALSE]
  row.names(observed) <- NULL
  best <- observed[1L, , drop = FALSE]

  best_settings <- .memwas_make_settings_for_pair(settings, settings$data,
                                                   best$L1_penalty[1L], best$L2_penalty[1L],
                                                   engine = backend,
                                                   autocor_regularization = autocor_regularization,
                                                   autocor_regularization_type = autocor_regularization_type,
                                                   autocor_regularization_alpha = autocor_regularization_alpha)

  final_fit <- NULL
  final_fit_error <- NULL
  if (isTRUE(refit_final)) {
    final_fit_try <- try(fit_MEMWAS(best_settings, verbose = FALSE, ...), silent = TRUE)
    if (inherits(final_fit_try, "try-error")) {
      final_fit_error <- as.character(final_fit_try)
      warning("Final MEMWAS refit failed: ", final_fit_error, call. = FALSE)
    } else {
      final_fit <- final_fit_try
    }
  } else {
    final_fit_error <- "Final model refit was skipped because `refit_final = FALSE`."
  }

  best_model <- .memwas_build_tuned_best_model(best = best,
                                               settings = best_settings,
                                               final_fit = final_fit,
                                               final_fit_error = final_fit_error)
  best <- .memwas_attach_best_model_columns(best, best_model)

  cv_details <- if (length(cv_detail_list)) do.call(rbind, cv_detail_list) else data.frame()
  row.names(cv_details) <- NULL

  out <- list(
    call = call,
    best = best,
    results = observed,
    cv_details = cv_details,
    best_fit = final_fit,
    best_model = best_model,
    best_settings = best_settings,
    final_fit_error = final_fit_error,
    fold_assignment = data.frame(row = seq_len(nrow(settings$data)), fold = row_folds,
                                 id = settings$data[[settings$id]], stringsAsFactors = FALSE),
    fold_assignment_ids = id_folds,
    settings = settings,
    tuning_control = list(
      K = K,
      approximation = settings$approximation,
      init_approximation = settings$init_approximation,
      se_method = settings$se_method,
      metric = metric,
      stability_metric = stability_metric,
      surrogate = surrogate,
      use_stability_penalty = use_stability_penalty,
      lambda_log_range = lambda_log_range,
      lambda_base = lambda_base,
      alpha_range = alpha_range,
      engine = backend,
      initial_design = initial_design,
      n_initial = n_initial,
      max_iter = max_iter,
      improvement_tol = improvement_tol,
      patience = patience,
      candidate_grid_size = candidate_grid_size,
      exploration = exploration,
      stability_weight = stability_weight,
      autocor_regularization = autocor_regularization,
      autocor_regularization_type = autocor_regularization_type,
      autocor_regularization_alpha = autocor_regularization_alpha
    ),
    convergence = list(completed_iterations = completed_iter,
                       total_evaluations = eval_counter,
                       stop_reason = stop_reason),
    fold_fits = if (isTRUE(keep_fold_fits)) fold_fit_list else NULL,
    note = paste(
      "Grouped K-fold cross-validation keeps all rows for each subject ID in the same fold.",
      "The final row in `best` minimizes median CV error plus the optional stability penalty."
    )
  )
  class(out) <- "tune_MEMWAS"
  out
}


# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# S3 METHOD: print summary for tune_MEMWAS objects -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
