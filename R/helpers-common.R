# Internal utility helpers shared by set_MEMWAS(), fit_MEMWAS(), and tune_MEMWAS().
# These helpers are intentionally not exported.

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

.memwas_supported_families <- function() {
  c("gaussian", "binomial", "poisson", "negative_binomial", "gamma", "exponential")
}

.memwas_normalize_family <- function(family) {
  if (inherits(family, "family")) family <- family$family
  family <- tolower(as.character(family)[1L])
  family <- gsub("\\s+", "_", family)
  family <- gsub("[-.]", "_", family)
  if (grepl("^negative_binomial", family)) family <- "negative_binomial"
  aliases <- c(normal = "gaussian", gaussian = "gaussian",
               logit = "binomial", logistic = "binomial", bernoulli = "binomial",
               binomial = "binomial",
               count = "poisson", poisson = "poisson",
               negative_binomial = "negative_binomial", negbin = "negative_binomial",
               neg_binomial = "negative_binomial", nb = "negative_binomial",
               nb2 = "negative_binomial", nbinom = "negative_binomial",
               gamma = "gamma", Gamma = "gamma",
               exponential = "exponential", exp = "exponential")
  if (!family %in% names(aliases)) {
    stop("Unsupported `family`. Supported values are 'gaussian', 'binomial', 'poisson', ",
         "'negative_binomial', 'gamma', and 'exponential'.",
         call. = FALSE)
  }
  unname(aliases[family])
}

.memwas_supported_approximations <- function() {
  c("laplace", "adaptive_gauss_hermite_quadrature", "adaptive_gaussian_quadrature",
    "saddlepoint", "skew_corrected_laplace", "variational_inference", "pql")
}

.memwas_normalize_approximation <- function(approximation = "laplace") {
  if (is.null(approximation) || length(approximation) == 0L || is.na(approximation[1L]) ||
      !nzchar(as.character(approximation[1L]))) {
    approximation <- "laplace"
  }
  key <- tolower(trimws(as.character(approximation[1L])))
  key <- gsub("[[:space:]-]+", "_", key)
  key <- gsub("\\.+", "_", key)
  aliases <- c(
    auto = "laplace", automatic = "laplace", default = "laplace",
    vi = "variational_inference", variational = "variational_inference",
    variational_inference = "variational_inference", mean_field = "variational_inference",
    mean_field_variational = "variational_inference", mfvi = "variational_inference",
    laplace = "laplace", laplace_approximation = "laplace",
    saddlepoint = "saddlepoint", saddlepoint_approximation = "saddlepoint",
    spa = "saddlepoint",
    skew_laplace = "skew_corrected_laplace",
    skew_corrected_laplace = "skew_corrected_laplace",
    skew_corrected_laplace_approximation = "skew_corrected_laplace",
    scl = "skew_corrected_laplace",
    adaptive_gaussian_quadrature = "adaptive_gaussian_quadrature",
    adaptive_gaussian = "adaptive_gaussian_quadrature", agq = "adaptive_gaussian_quadrature",
    adaptive_gauss_hermite_quadrature = "adaptive_gauss_hermite_quadrature",
    adaptive_gauss_hermite = "adaptive_gauss_hermite_quadrature",
    adaptive_gh = "adaptive_gauss_hermite_quadrature", aghq = "adaptive_gauss_hermite_quadrature",
    pql = "pql", penalized_quasi_likelihood = "pql"
  )
  if (!key %in% names(aliases)) {
    stop("Unsupported `approximation`. Supported values are 'laplace', ",
         "'adaptive_gauss_hermite_quadrature', 'adaptive_gaussian_quadrature', ",
         "'saddlepoint', 'skew_corrected_laplace', 'variational_inference', and legacy 'pql'.",
         call. = FALSE)
  }
  unname(aliases[key])
}

.memwas_approximation_label <- function(approximation) {
  if (!is.null(approximation) && length(approximation) > 0L &&
      tolower(as.character(approximation[1L])) %in% c("exact_gaussian", "gaussian_exact")) {
    return("Exact Gaussian marginal likelihood")
  }
  approximation <- .memwas_normalize_approximation(approximation)
  switch(approximation,
         variational_inference = "Variational inference",
         laplace = "Laplace approximation",
         saddlepoint = "Saddlepoint approximation",
         skew_corrected_laplace = "Skew-corrected Laplace approximation",
         adaptive_gaussian_quadrature = "Adaptive Gaussian quadrature",
         adaptive_gauss_hermite_quadrature = "Adaptive Gauss-Hermite quadrature",
         pql = "Penalized quasi-likelihood",
         approximation)
}

.memwas_parse_random_formula <- function(random) {
  if (is.null(random)) return(stats::as.formula("~ 1"))
  txt <- if (inherits(random, "formula")) paste(deparse(random), collapse = "") else as.character(random)[1L]
  txt <- trimws(txt)
  if (grepl("\\|", txt)) {
    left <- strsplit(txt, "\\|")[[1L]][1L]
    left <- gsub("^\\s*~", "", left)
    txt <- paste("~", left)
  }
  if (!grepl("^\\s*~", txt)) txt <- paste("~", txt)
  stats::as.formula(txt)
}

.memwas_add_terms_to_formula <- function(formula, add_terms) {
  add_terms <- unique(add_terms[nzchar(add_terms)])
  if (length(add_terms) == 0L) return(formula)
  trm <- stats::terms(formula)
  response <- deparse(formula[[2L]])
  labels <- attr(trm, "term.labels")
  intercept <- attr(trm, "intercept") == 1L
  new_labels <- unique(c(labels, add_terms))
  stats::reformulate(new_labels, response = response, intercept = intercept,
                     env = environment(formula))
}

.memwas_make_default_settings <- function(call, verbose = TRUE) {
  list(
    call = call,
    formula = NULL,
    formal_formula = NULL,
    family = "gaussian",
    data = NULL,
    id = NULL,
    time = NULL,
    random = stats::as.formula("~ 1"),
    random_cov = "diagonal",
    autocor = "AR(1)",
    serial = NULL,
    L1_penalty = 0,
    L2_penalty = 0,
    control = list(),
    method = "ML",
    approximation = "laplace",
    init_approximation = "variational_inference",
    se_method = "hessian",
    dot_predictors = NULL,
    dot_alternative = NULL,
    dot_threshold = 0,
    dot_alpha = 0.05,
    dot_spec = NULL,
    engine = "R",
    autocorrelation_check = "All",
    distribution_link_check = "All",
    conditional_independence_check = "All",
    random_effects_normality_check = "All",
    random_effects_predictor_independence_check = "All",
    homogeneity_variance_check = "All",
    assumption_check_spec = NULL,
    assumption_check_methods = data.frame(),
    assumption_checks = list(),
    assumption_check_model = NA_character_,
    nonlinear_summary = data.frame(),
    turning_points = list(),
    spline_variables = character(0L),
    spline_info = list(),
    all_screened_spline_info = list(),
    baseline_screen_metrics = NULL,
    extra = list(),
    verbose = verbose,
    note = paste(
      "Gaussian models use direct marginal ML/REML.",
      "Non-Gaussian models use variational-inference initialization followed by the selected native-C++ marginal approximation; Laplace is the default final approximation."
    )
  )
}

.memwas_check_scalar_nonnegative <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0) {
    stop("`", nm, "` must be a single non-negative numeric value.", call. = FALSE)
  }
}

.memwas_clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}


.memwas_dot_default_names <- function() c(".default", "default", ".all", "all", "*")

.memwas_is_dot_default_name <- function(x) {
  tolower(trimws(as.character(x))) %in% .memwas_dot_default_names()
}

.memwas_model_terms_from_assign <- function(coef_terms, assign = NULL, term_labels = NULL) {
  coef_terms <- as.character(coef_terms)
  out <- coef_terms
  if (is.null(assign) || length(assign) != length(coef_terms)) return(out)
  assign <- suppressWarnings(as.integer(assign))
  term_labels <- as.character(term_labels %||% character(0L))
  for (i in seq_along(out)) {
    ai <- assign[i]
    if (is.na(ai)) next
    if (ai == 0L) {
      out[i] <- "(Intercept)"
    } else if (ai >= 1L && ai <= length(term_labels)) {
      out[i] <- term_labels[ai]
    }
  }
  out
}

.memwas_dot_vector <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x)) return(unlist(x, use.names = TRUE))
  x
}

.memwas_normalize_dot_alternative <- function(x) {
  if (length(x) == 0L) return(character(0L))
  raw <- trimws(tolower(as.character(x)))
  key <- gsub("[[:space:]_-]+", ".", raw)
  key <- gsub("\\.+", ".", key)
  out <- rep(NA_character_, length(key))
  out[key %in% c("greater", "gt", "positive", "pos", "above", "right", ">", "larger")] <- "greater"
  out[key %in% c("less", "lt", "negative", "neg", "below", "left", "<", "smaller")] <- "less"
  out[key %in% c("two.sided", "twosided", "two.side", "two.sides", "two", "both", "two.tailed", "two.tail", "=", "!=")] <- "two.sided"
  if (anyNA(out)) {
    bad <- unique(raw[is.na(out)])
    stop("`dot_alternative` values must be 'greater', 'less', or 'two.sided'. Invalid value(s): ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}

.memwas_dot_miv <- function(dot_threshold, alternative = NULL) {
  threshold <- suppressWarnings(as.numeric(dot_threshold))
  if (any(!is.finite(threshold))) stop("`dot_threshold` must contain finite numeric value(s).", call. = FALSE)
  if (any(threshold < 0)) {
    warning("`dot_threshold` is interpreted as a non-negative minimally important value; negative value(s) were converted to absolute value(s).", call. = FALSE)
  }
  abs(threshold)
}

.memwas_dot_null_value <- function(dot_threshold, alternative) {
  threshold <- .memwas_dot_miv(dot_threshold, alternative)
  alternative <- .memwas_normalize_dot_alternative(alternative)
  out <- rep(0, length(threshold))
  out[alternative == "greater"] <- threshold[alternative == "greater"]
  out[alternative == "less"] <- -threshold[alternative == "less"]
  # For two-sided minimum-effect tests the null region is [-threshold, +threshold],
  # so there is no single signed boundary. Keep null_value at 0 and use
  # threshold/minimal_important_value for the p-value calculation.
  as.numeric(out)
}

.memwas_dot_hypothesis <- function(alternative, threshold) {
  alternative <- .memwas_normalize_dot_alternative(alternative)
  threshold <- .memwas_dot_miv(threshold, alternative)
  out <- character(length(alternative))
  greater <- alternative == "greater"
  less <- alternative == "less"
  two <- alternative == "two.sided"
  out[greater] <- paste0("H1: beta > ", format(signif(threshold[greater], 6L), trim = TRUE, scientific = FALSE))
  out[less & threshold == 0] <- "H1: beta < 0"
  out[less & threshold > 0] <- paste0("H1: beta < -", format(signif(threshold[less & threshold > 0], 6L), trim = TRUE, scientific = FALSE))
  out[two & threshold == 0] <- "H1: beta != 0"
  out[two & threshold > 0] <- paste0("H1: |beta| > ", format(signif(threshold[two & threshold > 0], 6L), trim = TRUE, scientific = FALSE))
  out
}

.memwas_dot_dataframe_to_components <- function(dot_alternative, dot_threshold, dot_alpha) {
  dfs <- list()
  if (is.data.frame(dot_alternative)) dfs$alternative <- dot_alternative
  if (is.data.frame(dot_threshold)) dfs$threshold <- dot_threshold
  if (is.data.frame(dot_alpha)) dfs$alpha <- dot_alpha

  if (length(dfs) == 0L) {
    return(list(dot_predictors = NULL, dot_alternative = dot_alternative,
                dot_threshold = dot_threshold, dot_alpha = dot_alpha))
  }

  predictors <- character(0L)
  alt <- NULL
  thr <- NULL
  alp <- NULL

  for (src in names(dfs)) {
    df <- dfs[[src]]
    nms <- names(df)
    term_col <- intersect(c("term", "predictor", "coefficient", "variable"), nms)[1L]
    if (is.na(term_col)) {
      stop("A DOT data frame must contain a `term` or `predictor` column.", call. = FALSE)
    }
    these_predictors <- as.character(df[[term_col]])
    if (anyNA(these_predictors) || any(!nzchar(these_predictors))) {
      stop("DOT predictor names must be non-missing, non-empty character values.", call. = FALSE)
    }
    predictors <- unique(c(predictors, these_predictors))

    alt_col <- intersect(c("alternative", "dot_alternative", "direction"), nms)[1L]
    threshold_col <- intersect(c("dot_threshold", "threshold", "minimal_important_value",
                                 "minimally_important_value", "minimal_value", "miv", "margin"), nms)[1L]
    alpha_col <- intersect(c("dot_alpha", "alpha", "significance_alpha", "significance_level"), nms)[1L]

    if (src == "alternative" && is.na(alt_col)) {
      stop("A data-frame `dot_alternative` must contain an `alternative` column.", call. = FALSE)
    }
    if (src == "threshold" && is.na(threshold_col)) {
      stop("A data-frame `dot_threshold` must contain a `threshold`, `dot_threshold`, or `minimal_important_value` column.", call. = FALSE)
    }
    if (src == "alpha" && is.na(alpha_col)) {
      stop("A data-frame `dot_alpha` must contain an `alpha` or `dot_alpha` column.", call. = FALSE)
    }

    if (!is.na(alt_col)) alt <- c(alt, setNames(as.character(df[[alt_col]]), these_predictors))
    if (!is.na(threshold_col)) thr <- c(thr, setNames(as.numeric(df[[threshold_col]]), these_predictors))
    if (!is.na(alpha_col)) alp <- c(alp, setNames(as.numeric(df[[alpha_col]]), these_predictors))
  }

  list(
    dot_predictors = predictors,
    dot_alternative = if (!is.null(alt)) alt else if (is.data.frame(dot_alternative)) NULL else dot_alternative,
    dot_threshold = if (!is.null(thr)) thr else if (is.data.frame(dot_threshold)) NULL else dot_threshold,
    dot_alpha = if (!is.null(alp)) alp else if (is.data.frame(dot_alpha)) NULL else dot_alpha
  )
}

.memwas_dot_component_parts <- function(x, default, label,
                                        cast = c("character", "numeric"),
                                        allow_default = TRUE) {
  cast <- match.arg(cast)
  make_value <- function(z) {
    if (cast == "numeric") {
      out <- suppressWarnings(as.numeric(z))
      if (anyNA(out)) stop("`", label, "` must contain numeric value(s).", call. = FALSE)
      out
    } else {
      as.character(z)
    }
  }

  out <- list(default = default, has_default = FALSE,
              named = if (cast == "numeric") numeric(0L) else character(0L),
              unnamed = if (cast == "numeric") numeric(0L) else character(0L))
  if (is.null(x)) return(out)
  x <- .memwas_dot_vector(x)
  if (is.data.frame(x)) stop("Internal DOT parsing error: unexpected data frame.", call. = FALSE)
  if (length(x) == 0L) return(out)

  nms <- names(x)
  if (is.null(nms)) nms <- rep("", length(x))
  nms <- as.character(nms)
  default_key <- nzchar(nms) & .memwas_is_dot_default_name(nms)
  if (any(default_key)) {
    if (sum(default_key) > 1L) stop("`", label, "` may contain only one .default value.", call. = FALSE)
    out$default <- make_value(x[default_key])[1L]
    out$has_default <- TRUE
  }

  named <- nzchar(nms) & !default_key
  if (any(named)) {
    vals <- make_value(x[named])
    names(vals) <- nms[named]
    out$named <- vals
  }

  unnamed <- !nzchar(nms) & !default_key
  if (any(unnamed)) {
    vals <- make_value(x[unnamed])
    if (length(vals) == 1L && isTRUE(allow_default)) {
      out$default <- vals[1L]
      out$has_default <- TRUE
    } else {
      out$unnamed <- vals
    }
  }
  out
}

.memwas_dot_unique_terms <- function(dot_predictors = NULL, ...) {
  terms <- character(0L)
  if (!is.null(dot_predictors)) terms <- c(terms, as.character(dot_predictors))
  parts <- list(...)
  for (part in parts) {
    if (is.null(part) || is.null(part$named)) next
    nms <- names(part$named)
    if (!is.null(nms)) terms <- c(terms, nms[nzchar(nms) & !.memwas_is_dot_default_name(nms)])
  }
  unique(terms[nzchar(terms) & !is.na(terms) & !.memwas_is_dot_default_name(terms)])
}

.memwas_dot_values_for_terms <- function(part, terms, explicit_predictors = NULL,
                                         default, label, normalize = NULL) {
  values <- rep(default, length(terms))
  names(values) <- terms
  if (!is.null(part) && !is.null(part$named) && length(part$named) > 0L) {
    nms <- names(part$named)
    keep <- nms %in% terms
    if (any(keep)) values[nms[keep]] <- part$named[keep]
  }
  if (!is.null(part) && !is.null(part$unnamed) && length(part$unnamed) > 0L) {
    if (is.null(explicit_predictors) || length(explicit_predictors) == 0L) {
      stop("Unnamed `", label, "` values longer than one require `dot_predictors` so they can be aligned by position.", call. = FALSE)
    }
    explicit_predictors <- as.character(explicit_predictors)
    if (length(part$unnamed) != length(explicit_predictors)) {
      stop("Unnamed `", label, "` values must have length 1, be named by predictor, or match the length of `dot_predictors`.", call. = FALSE)
    }
    idx <- match(explicit_predictors, terms)
    ok <- !is.na(idx)
    if (any(ok)) values[idx[ok]] <- part$unnamed[ok]
  }
  if (!is.null(normalize)) values <- normalize(values)
  values
}

.memwas_validate_dot_settings <- function(dot_predictors = NULL,
                                          dot_alternative = NULL,
                                          dot_threshold = 0,
                                          dot_alpha = 0.05) {
  df_parts <- .memwas_dot_dataframe_to_components(dot_alternative, dot_threshold, dot_alpha)
  if (!is.null(df_parts$dot_predictors)) dot_predictors <- df_parts$dot_predictors
  dot_alternative <- df_parts$dot_alternative
  dot_threshold <- df_parts$dot_threshold
  dot_alpha <- df_parts$dot_alpha

  alt_part <- .memwas_dot_component_parts(dot_alternative, "two.sided", "dot_alternative",
                                          cast = "character", allow_default = TRUE)
  threshold_part <- .memwas_dot_component_parts(dot_threshold, 0, "dot_threshold",
                                                cast = "numeric", allow_default = TRUE)
  alpha_part <- .memwas_dot_component_parts(dot_alpha, 0.05, "dot_alpha",
                                            cast = "numeric", allow_default = TRUE)

  default_alternative <- .memwas_normalize_dot_alternative(alt_part$default)[1L]
  default_dot_threshold <- .memwas_dot_miv(as.numeric(threshold_part$default)[1L], default_alternative)
  default_alpha <- as.numeric(alpha_part$default)[1L]
  if (!is.finite(default_alpha) || default_alpha <= 0 || default_alpha >= 1) {
    stop("`dot_alpha` must contain numeric value(s) strictly between 0 and 1.", call. = FALSE)
  }

  default_threshold <- .memwas_dot_miv(default_dot_threshold, default_alternative)[1L]
  default_null_value <- .memwas_dot_null_value(default_dot_threshold, default_alternative)[1L]
  predictors <- .memwas_dot_unique_terms(dot_predictors, alt_part, threshold_part, alpha_part)
  if (length(predictors) == 0L &&
      (length(alt_part$unnamed) > 0L || length(threshold_part$unnamed) > 0L || length(alpha_part$unnamed) > 0L)) {
    stop("Unnamed DOT vectors with length greater than 1 require `dot_predictors` or names on the vector values.", call. = FALSE)
  }

  if (length(predictors) == 0L) {
    tab <- data.frame(term = character(0L), alternative = character(0L),
                      dot_threshold = numeric(0L), threshold = numeric(0L),
                      minimal_important_value = numeric(0L),
                      null_value = numeric(0L), test_threshold = numeric(0L),
                      null_boundary = numeric(0L), alpha = numeric(0L), hypothesis = character(0L),
                      stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    alt <- .memwas_dot_values_for_terms(alt_part, predictors,
                                        explicit_predictors = dot_predictors,
                                        default = default_alternative,
                                        label = "dot_alternative",
                                        normalize = .memwas_normalize_dot_alternative)
    dot_thr <- .memwas_dot_values_for_terms(threshold_part, predictors,
                                            explicit_predictors = dot_predictors,
                                            default = default_dot_threshold,
                                            label = "dot_threshold")
    alp <- .memwas_dot_values_for_terms(alpha_part, predictors,
                                        explicit_predictors = dot_predictors,
                                        default = default_alpha,
                                        label = "dot_alpha")
    dot_thr <- .memwas_dot_miv(dot_thr, alt)
    if (any(!is.finite(alp) | alp <= 0 | alp >= 1)) {
      stop("`dot_alpha` must contain numeric value(s) strictly between 0 and 1.", call. = FALSE)
    }
    threshold <- .memwas_dot_miv(dot_thr, alt)
    null_value <- .memwas_dot_null_value(dot_thr, alt)
    tab <- data.frame(term = predictors, alternative = as.character(alt),
                      dot_threshold = as.numeric(threshold),
                      threshold = as.numeric(threshold),
                      minimal_important_value = as.numeric(threshold),
                      null_value = as.numeric(null_value),
                      test_threshold = as.numeric(null_value),
                      null_boundary = as.numeric(null_value),
                      alpha = as.numeric(alp),
                      hypothesis = .memwas_dot_hypothesis(alt, threshold),
                      stringsAsFactors = FALSE, check.names = FALSE)
    row.names(tab) <- NULL
  }

  out <- list(
    predictors = predictors,
    table = tab,
    default_alternative = as.character(default_alternative)[1L],
    default_dot_threshold = as.numeric(default_dot_threshold)[1L],
    default_threshold = as.numeric(default_threshold)[1L],
    default_minimal_important_value = as.numeric(default_threshold)[1L],
    default_null_value = as.numeric(default_null_value)[1L],
    default_test_threshold = as.numeric(default_null_value)[1L],
    default_null_boundary = as.numeric(default_null_value)[1L],
    default_alpha = as.numeric(default_alpha)[1L],
    default_hypothesis = .memwas_dot_hypothesis(default_alternative, default_threshold),
    threshold_scale = "minimal_important_value",
    note = paste(
      "`dot_threshold` is interpreted as a non-negative minimally important value for each predictor.",
      "For alternative = 'greater', MEMWAS tests H1: beta > +dot_threshold; for alternative = 'less', MEMWAS tests H1: beta < -dot_threshold.",
      "For alternative = 'two.sided', dot_threshold = 0 gives the usual two-sided Wald test and dot_threshold > 0 gives a two-sided minimum-effect Wald test for |beta| > dot_threshold.",
      "Summary tables report `threshold`/`minimal_important_value` as the minimally important value and `null_value`/`test_threshold`/`null_boundary` as the signed one-tailed boundary used for greater/less tests."
    )
  )
  class(out) <- "MEMWAS_dot"
  out
}

.memwas_dot_term_metadata <- function(terms, term_map = NULL, model_matrix = NULL, formula = NULL) {
  terms <- as.character(terms)
  model_term <- terms
  if (is.data.frame(term_map) && all(c("term", "model_term") %in% names(term_map))) {
    idx <- match(terms, as.character(term_map$term))
    ok <- !is.na(idx)
    if (any(ok)) model_term[ok] <- as.character(term_map$model_term)[idx[ok]]
  } else {
    assign <- if (!is.null(model_matrix)) attr(model_matrix, "assign") else NULL
    labels <- character(0L)
    if (inherits(formula, "formula")) {
      labels <- attr(stats::terms(formula), "term.labels") %||% character(0L)
    }
    model_term <- .memwas_model_terms_from_assign(terms, assign = assign, term_labels = labels)
  }
  data.frame(term = terms, model_term = as.character(model_term),
             stringsAsFactors = FALSE, check.names = FALSE)
}

.memwas_dot_config_for_terms <- function(terms, dot_spec = NULL, term_map = NULL,
                                         model_matrix = NULL, formula = NULL) {
  if (is.null(dot_spec)) dot_spec <- .memwas_validate_dot_settings()
  meta <- .memwas_dot_term_metadata(terms, term_map = term_map,
                                    model_matrix = model_matrix, formula = formula)
  terms <- meta$term

  alternative <- rep(dot_spec$default_alternative %||% "two.sided", length(terms))
  threshold <- rep(dot_spec$default_threshold %||% 0, length(terms))
  null_value <- rep(dot_spec$default_null_value %||% 0, length(terms))
  alpha <- rep(dot_spec$default_alpha %||% 0.05, length(terms))
  matched_by <- rep("default", length(terms))

  is_intercept <- terms %in% c("(Intercept)", "Intercept") | meta$model_term %in% c("(Intercept)", "Intercept")
  if (any(is_intercept)) {
    alternative[is_intercept] <- "two.sided"
    threshold[is_intercept] <- 0
    null_value[is_intercept] <- 0
    alpha[is_intercept] <- 0.05
    matched_by[is_intercept] <- "intercept_default"
  }

  if (!is.null(dot_spec$table) && nrow(dot_spec$table) > 0L) {
    table_test_threshold <- function(i) {
      if ("test_threshold" %in% names(dot_spec$table)) return(dot_spec$table$test_threshold[i])
      if ("null_value" %in% names(dot_spec$table)) return(dot_spec$table$null_value[i])
      .memwas_dot_null_value(dot_spec$table$threshold[i], dot_spec$table$alternative[i])
    }
    for (i in seq_len(nrow(dot_spec$table))) {
      idx <- which(meta$model_term == dot_spec$table$term[i])
      if (length(idx) > 0L) {
        alternative[idx] <- dot_spec$table$alternative[i]
        threshold[idx] <- dot_spec$table$threshold[i]
        null_value[idx] <- table_test_threshold(i)
        alpha[idx] <- dot_spec$table$alpha[i]
        matched_by[idx] <- "model_term"
      }
    }
    for (i in seq_len(nrow(dot_spec$table))) {
      idx <- which(terms == dot_spec$table$term[i])
      if (length(idx) > 0L) {
        alternative[idx] <- dot_spec$table$alternative[i]
        threshold[idx] <- dot_spec$table$threshold[i]
        null_value[idx] <- table_test_threshold(i)
        alpha[idx] <- dot_spec$table$alpha[i]
        matched_by[idx] <- "coefficient"
      }
    }
  }

  data.frame(term = terms, model_term = meta$model_term,
             dot_target = ifelse(matched_by == "coefficient", terms, meta$model_term),
             alternative = alternative,
             threshold = as.numeric(threshold),
             minimal_important_value = as.numeric(threshold),
             null_value = as.numeric(null_value),
             test_threshold = as.numeric(null_value),
             null_boundary = as.numeric(null_value),
             alpha = as.numeric(alpha),
             matched_by = matched_by,
             hypothesis = .memwas_dot_hypothesis(alternative, threshold),
             stringsAsFactors = FALSE, check.names = FALSE)
}

.memwas_apply_dot_to_coefficient_table <- function(coef_table, dot_spec = NULL,
                                                   term_map = NULL, model_matrix = NULL,
                                                   formula = NULL) {
  if (is.null(coef_table) || !is.data.frame(coef_table) || nrow(coef_table) == 0L) return(coef_table)
  if (!"term" %in% names(coef_table)) stop("Coefficient table must contain a `term` column.", call. = FALSE)
  if (!"estimate" %in% names(coef_table)) stop("Coefficient table must contain an `estimate` column.", call. = FALSE)
  if (!"std_error" %in% names(coef_table)) stop("Coefficient table must contain a `std_error` column.", call. = FALSE)

  if (is.null(term_map) && "model_term" %in% names(coef_table)) {
    term_map <- data.frame(term = coef_table$term, model_term = coef_table$model_term,
                           stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.null(term_map) && "predictor" %in% names(coef_table)) {
    term_map <- data.frame(term = coef_table$term, model_term = coef_table$predictor,
                           stringsAsFactors = FALSE, check.names = FALSE)
  }

  cfg <- .memwas_dot_config_for_terms(coef_table$term, dot_spec, term_map = term_map,
                                      model_matrix = model_matrix, formula = formula)
  estimate <- as.numeric(coef_table$estimate)
  std_error <- as.numeric(coef_table$std_error)
  threshold <- as.numeric(cfg$threshold)
  null_value <- as.numeric(cfg$null_value)
  statistic <- rep(NA_real_, length(estimate))
  ok_base <- is.finite(estimate) & is.finite(std_error) & std_error > 0

  p_value <- rep(NA_real_, length(statistic))
  greater <- ok_base & is.finite(null_value) & cfg$alternative == "greater"
  less <- ok_base & is.finite(null_value) & cfg$alternative == "less"
  two <- ok_base & is.finite(threshold) & cfg$alternative == "two.sided"

  statistic[greater] <- (estimate[greater] - null_value[greater]) / std_error[greater]
  statistic[less] <- (estimate[less] - null_value[less]) / std_error[less]
  p_value[greater] <- stats::pnorm(statistic[greater], lower.tail = FALSE)
  p_value[less] <- stats::pnorm(statistic[less], lower.tail = TRUE)

  statistic[two] <- ifelse(threshold[two] == 0,
                           estimate[two] / std_error[two],
                           (abs(estimate[two]) - threshold[two]) / std_error[two])
  p_value[two] <- ifelse(threshold[two] == 0,
                         2 * stats::pnorm(abs(estimate[two] / std_error[two]), lower.tail = FALSE),
                         pmin(1, 2 * stats::pnorm(statistic[two], lower.tail = FALSE)))
  significant <- ifelse(is.finite(p_value) & p_value <= cfg$alpha, TRUE, FALSE)

  coef_table$model_term <- cfg$model_term
  coef_table$dot_target <- cfg$dot_target
  coef_table$dot_threshold <- threshold
  coef_table$threshold <- threshold
  coef_table$minimal_important_value <- cfg$minimal_important_value
  coef_table$null_value <- null_value
  coef_table$test_threshold <- null_value
  coef_table$null_boundary <- null_value
  coef_table$alternative <- cfg$alternative
  coef_table$direction <- ifelse(cfg$alternative == "greater", "positive",
                                 ifelse(cfg$alternative == "less", "negative", "two-sided"))
  coef_table$alpha <- cfg$alpha
  coef_table$hypothesis <- cfg$hypothesis
  coef_table$test <- ifelse(cfg$alternative == "two.sided" & threshold > 0,
                            "two-sided minimum-effect Wald",
                            ifelse(cfg$alternative == "two.sided", "two-sided Wald",
                                   "directional one-tailed Wald vs minimally important value"))
  coef_table$statistic <- as.numeric(statistic)
  coef_table$p_value <- as.numeric(p_value)
  coef_table$significant <- as.logical(significant)
  coef_table$DOT_label <- ifelse(cfg$alternative == "two.sided", "two-sided", "DOT")
  coef_table$dot_matched_by <- cfg$matched_by

  desired <- c("term", "model_term", "dot_target", "estimate", "std_error",
               "dot_threshold", "threshold", "minimal_important_value", "null_value", "test_threshold", "null_boundary",
               "alternative", "direction", "alpha", "hypothesis", "test",
               "statistic", "p_value", "significant", "DOT_label", "dot_matched_by")
  coef_table[, c(intersect(desired, names(coef_table)), setdiff(names(coef_table), desired)), drop = FALSE]
}

.memwas_warn_unmatched_dot_predictors <- function(dot_spec, coefficient_terms,
                                                 predictor_terms = NULL,
                                                 term_map = NULL) {
  if (is.null(dot_spec) || is.null(dot_spec$predictors) || length(dot_spec$predictors) == 0L) {
    return(invisible(character(0L)))
  }
  coefficient_terms <- as.character(coefficient_terms)
  if (is.data.frame(term_map) && "model_term" %in% names(term_map)) {
    predictor_terms <- unique(c(predictor_terms, as.character(term_map$model_term)))
  }
  predictor_terms <- if (is.null(predictor_terms)) character(0L) else as.character(predictor_terms)
  matched <- dot_spec$predictors %in% coefficient_terms | dot_spec$predictors %in% predictor_terms
  unmatched <- dot_spec$predictors[!matched]
  if (length(unmatched) > 0L) {
    warning("DOT predictor(s) were not matched to fitted fixed-effect coefficient terms or model terms: ",
            paste(unmatched, collapse = ", "),
            ". Use exact coefficient names from the fixed-effect summary table or predictor/model-term names from the model formula.",
            call. = FALSE)
  }
  invisible(unmatched)
}
