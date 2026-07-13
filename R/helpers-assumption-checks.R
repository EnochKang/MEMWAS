# Internal helpers for base-R GLMM assumption screening.
# These functions are intentionally not exported and do not add package dependencies.

.memwas_assumption_method_catalog <- function() {
  list(
    autocorrelation_check = c("DurbinWatson", "LjungBox", "Lag1Correlation", "Runs"),
    distribution_link_check = c("PearsonDispersion", "DevianceDispersion", "LinkTest", "QuantileResidualNormality", "GroupedCalibration"),
    conditional_independence_check = c("WithinClusterLag1", "ClusterMeanResidual", "Runs"),
    random_effects_normality_check = c("Shapiro", "JarqueBera", "SkewKurtosis"),
    random_effects_predictor_independence_check = c("GroupMeanAssociation", "Correlation", "RankCorrelation"),
    homogeneity_variance_check = c("BreuschPagan", "White", "LeveneFitted", "LeveneGroup")
  )
}

.memwas_assumption_method_key <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

.memwas_assumption_alias_vector <- function(choices, extra = list()) {
  out <- stats::setNames(choices, .memwas_assumption_method_key(choices))
  if (length(extra) > 0L) {
    for (nm in names(extra)) {
      out[.memwas_assumption_method_key(extra[[nm]])] <- nm
    }
  }
  out
}

.memwas_assumption_alias_catalog <- function() {
  catalog <- .memwas_assumption_method_catalog()
  list(
    autocorrelation_check = .memwas_assumption_alias_vector(
      catalog$autocorrelation_check,
      list(DurbinWatson = c("DW", "Durbin-Watson", "Durbin Watson"),
           LjungBox = c("Ljung-Box", "Box-Ljung", "Portmanteau", "BoxPierce", "Box-Pierce"),
           Lag1Correlation = c("Lag1", "Lag-1", "Lag 1", "Lag1Cor", "ACF1"),
           Runs = c("RunsTest", "Runs Test", "WaldWolfowitz", "Wald-Wolfowitz"))
    ),
    distribution_link_check = .memwas_assumption_alias_vector(
      catalog$distribution_link_check,
      list(PearsonDispersion = c("Pearson", "Pearson Chi-square", "Overdispersion"),
           DevianceDispersion = c("Deviance", "Deviance Chi-square"),
           LinkTest = c("Link", "Pregibon", "Pregibon Link", "Link Function"),
           QuantileResidualNormality = c("QuantileResidual", "Quantile Residual", "DunnSmyth", "Dunn-Smyth"),
           GroupedCalibration = c("Calibration", "HosmerLemeshow", "Hosmer-Lemeshow", "Grouped Calibration"))
    ),
    conditional_independence_check = .memwas_assumption_alias_vector(
      catalog$conditional_independence_check,
      list(WithinClusterLag1 = c("WithinCluster", "Within Cluster Lag1", "Lag1", "Lag-1"),
           ClusterMeanResidual = c("ClusterMean", "Cluster Mean", "ResidualByCluster", "ANOVA"),
           Runs = c("RunsTest", "Runs Test", "WaldWolfowitz", "Wald-Wolfowitz"))
    ),
    random_effects_normality_check = .memwas_assumption_alias_vector(
      catalog$random_effects_normality_check,
      list(Shapiro = c("ShapiroWilk", "Shapiro-Wilk"),
           JarqueBera = c("JB", "Jarque-Bera"),
           SkewKurtosis = c("Skew", "Kurtosis", "SkewnessKurtosis", "Skewness-Kurtosis"))
    ),
    random_effects_predictor_independence_check = .memwas_assumption_alias_vector(
      catalog$random_effects_predictor_independence_check,
      list(GroupMeanAssociation = c("Mundlak", "CRE", "CorrelatedRandomEffects", "Group Means"),
           Correlation = c("Pearson", "PearsonCorrelation", "Pearson Correlation"),
           RankCorrelation = c("Spearman", "SpearmanCorrelation", "Rank Correlation"))
    ),
    homogeneity_variance_check = .memwas_assumption_alias_vector(
      catalog$homogeneity_variance_check,
      list(BreuschPagan = c("BP", "Breusch-Pagan", "CookWeisberg", "Cook-Weisberg"),
           White = c("WhiteTest", "White Test"),
           LeveneFitted = c("LeveneByFitted", "FittedLevene", "BrownForsytheFitted"),
           LeveneGroup = c("LeveneByGroup", "GroupLevene", "BrownForsythe", "Brown-Forsythe"))
    )
  )
}

.memwas_assumption_normalize_methods <- function(x, arg, choices, aliases) {
  if (is.null(x)) return(character(0L))
  if (is.logical(x)) {
    if (length(x) != 1L || is.na(x)) stop("`", arg, "` must be TRUE, FALSE, 'All', 'None', or a supported method vector.", call. = FALSE)
    return(if (isTRUE(x)) choices else character(0L))
  }
  raw <- trimws(as.character(x))
  raw <- raw[nzchar(raw)]
  if (!length(raw)) return(character(0L))
  key <- .memwas_assumption_method_key(raw)
  all_keys <- c("all", "any", "everything", "full", "complete", "star")
  none_keys <- c("none", "no", "false", "off", "skip", "omit", "disable", "disabled")
  key[key == ""] <- "none"
  if (any(key %in% all_keys) || any(raw == "*")) return(choices)
  if (any(key %in% none_keys)) {
    if (length(key) == 1L) return(character(0L))
    keep <- !(key %in% none_keys)
    raw <- raw[keep]
    key <- key[keep]
  }
  out <- unname(aliases[key])
  bad <- raw[is.na(out)]
  if (length(bad) > 0L) {
    stop("Unsupported method(s) in `", arg, "`: ", paste(unique(bad), collapse = ", "),
         ". Supported values are 'All', 'None', or: ", paste(choices, collapse = ", "), ".",
         call. = FALSE)
  }
  unique(out)
}

.memwas_validate_assumption_check_settings <- function(
    autocorrelation_check = "All",
    distribution_link_check = "All",
    conditional_independence_check = "All",
    random_effects_normality_check = "All",
    random_effects_predictor_independence_check = "All",
    homogeneity_variance_check = "All") {
  catalog <- .memwas_assumption_method_catalog()
  aliases <- .memwas_assumption_alias_catalog()
  out <- list(
    autocorrelation_check = .memwas_assumption_normalize_methods(autocorrelation_check, "autocorrelation_check", catalog$autocorrelation_check, aliases$autocorrelation_check),
    distribution_link_check = .memwas_assumption_normalize_methods(distribution_link_check, "distribution_link_check", catalog$distribution_link_check, aliases$distribution_link_check),
    conditional_independence_check = .memwas_assumption_normalize_methods(conditional_independence_check, "conditional_independence_check", catalog$conditional_independence_check, aliases$conditional_independence_check),
    random_effects_normality_check = .memwas_assumption_normalize_methods(random_effects_normality_check, "random_effects_normality_check", catalog$random_effects_normality_check, aliases$random_effects_normality_check),
    random_effects_predictor_independence_check = .memwas_assumption_normalize_methods(random_effects_predictor_independence_check, "random_effects_predictor_independence_check", catalog$random_effects_predictor_independence_check, aliases$random_effects_predictor_independence_check),
    homogeneity_variance_check = .memwas_assumption_normalize_methods(homogeneity_variance_check, "homogeneity_variance_check", catalog$homogeneity_variance_check, aliases$homogeneity_variance_check)
  )
  class(out) <- "MEMWAS_assumption_spec"
  out
}

.memwas_assumption_is_empty_spec <- function(spec) {
  if (is.null(spec) || !is.list(spec)) return(TRUE)
  all(vapply(spec, function(z) length(z) == 0L, logical(1L)))
}

.memwas_assumption_decision <- function(p_value, alpha) {
  if (!is.finite(p_value)) return("Not tested")
  if (p_value <= alpha) "Potential violation" else "No evidence detected"
}

.memwas_assumption_result <- function(assumption, method, test, target = NA_character_,
                                      statistic = NA_real_, df = NA_real_, auxiliary = NA_real_,
                                      p_value = NA_real_, n = NA_integer_, alpha = 0.05,
                                      note = "") {
  data.frame(
    assumption = as.character(assumption)[1L],
    method = as.character(method)[1L],
    test = as.character(test)[1L],
    target = as.character(target)[1L],
    statistic = as.numeric(statistic)[1L],
    df = as.numeric(df)[1L],
    auxiliary = as.numeric(auxiliary)[1L],
    p_value = as.numeric(p_value)[1L],
    alpha = as.numeric(alpha)[1L],
    decision = .memwas_assumption_decision(as.numeric(p_value)[1L], as.numeric(alpha)[1L]),
    n = as.integer(n)[1L],
    note = as.character(note)[1L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.memwas_assumption_empty_result <- function(assumption, method, note, alpha = 0.05,
                                            test = method, target = NA_character_) {
  .memwas_assumption_result(assumption = assumption, method = method, test = test,
                            target = target, alpha = alpha, note = note)
}

.memwas_assumption_empty_table <- function() {
  data.frame(assumption = character(0L), method = character(0L), test = character(0L),
             target = character(0L), statistic = numeric(0L), df = numeric(0L),
             auxiliary = numeric(0L), p_value = numeric(0L), alpha = numeric(0L),
             decision = character(0L), n = integer(0L), note = character(0L),
             stringsAsFactors = FALSE, check.names = FALSE)
}

.memwas_split_assumption_table <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L || !"assumption" %in% names(tab)) return(list())
  split(tab, tab$assumption, drop = TRUE)
}

.memwas_assumption_alpha <- function(settings) {
  control <- settings$control %||% list()
  alpha <- suppressWarnings(as.numeric(control$assumption_alpha %||% control$diagnostic_alpha %||% 0.05)[1L])
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  alpha
}

.memwas_group_time_from_fit <- function(fit, settings, n) {
  d <- fit$data %||% settings$data %||% NULL
  group <- NULL
  time <- NULL
  if (is.data.frame(d) && nrow(d) == n) {
    if (!is.null(settings$id) && settings$id %in% names(d)) group <- d[[settings$id]]
    if (!is.null(settings$time) && settings$time %in% names(d)) time <- d[[settings$time]]
  }
  if (is.null(group) && is.list(fit$groups)) {
    group <- rep(NA_character_, n)
    for (g in seq_along(fit$groups)) group[fit$groups[[g]]] <- names(fit$groups)[g] %||% as.character(g)
  }
  if (is.null(time) && is.list(fit$groups)) {
    time <- rep(NA_real_, n)
    for (g in seq_along(fit$groups)) time[fit$groups[[g]]] <- seq_along(fit$groups[[g]])
  }
  list(group = group, time = time)
}

.memwas_assumption_family_variance <- function(family, mu, fit, settings) {
  family <- .memwas_normalize_family(family)
  mu <- as.numeric(mu)
  if (family == "gaussian") {
    sig <- suppressWarnings(as.numeric(fit$residual_sigma %||% NA_real_)[1L])
    if (!is.finite(sig) || sig <= 0) sig <- stats::sd(as.numeric(fit$residuals %||% rep(NA_real_, length(mu))), na.rm = TRUE)
    if (!is.finite(sig) || sig <= 0) sig <- 1
    return(rep(sig^2, length(mu)))
  }
  if (family == "binomial") return(pmax(mu * (1 - mu), 1e-8))
  if (family == "poisson") return(pmax(mu, 1e-8))
  if (family == "negative_binomial") {
    theta <- suppressWarnings(as.numeric(fit$family_parameters$theta %||% fit$family_parameters$size %||% settings$control$negative_binomial_theta %||% settings$control$nb_theta %||% NA_real_)[1L])
    if (!is.finite(theta) || theta <= 0) theta <- .memwas_estimate_nb_theta(fit$y %||% NULL, settings$control %||% list())
    return(pmax(mu + mu^2 / pmax(theta, 1e-8), 1e-8))
  }
  if (family == "gamma") {
    shape <- suppressWarnings(as.numeric(fit$family_parameters$shape %||% settings$control$gamma_shape %||% NA_real_)[1L])
    if (!is.finite(shape) || shape <= 0) shape <- .memwas_estimate_gamma_shape(fit$y %||% NULL, settings$control %||% list())
    return(pmax(mu^2 / pmax(shape, 1e-8), 1e-8))
  }
  if (family == "exponential") return(pmax(mu^2, 1e-8))
  rep(1, length(mu))
}

.memwas_conditional_fitted <- function(fit, settings) {
  y <- as.numeric(fit$y %||% numeric(0L))
  X <- fit$X %||% matrix(numeric(0L), nrow = length(y), ncol = 0L)
  Z <- fit$Z %||% matrix(numeric(0L), nrow = length(y), ncol = 0L)
  beta <- as.numeric(fit$coefficients %||% numeric(0L))
  if (!length(y) || !length(beta) || ncol(X) != length(beta)) {
    return(list(y = y, fitted_link = as.numeric(fit$fitted_link %||% fit$fitted %||% rep(NA_real_, length(y))),
                fitted_response = as.numeric(fit$fitted_response %||% fit$fitted %||% rep(NA_real_, length(y)))))
  }
  fixed_link <- as.vector(X %*% beta)
  random_link <- rep(0, length(y))
  if (ncol(Z) > 0L && is.list(fit$groups) && length(fit$groups) > 0L && is.list(fit$random_effects)) {
    for (g in seq_along(fit$groups)) {
      ii <- fit$groups[[g]]
      b <- fit$random_effects[[g]]
      if (is.null(b) && !is.null(names(fit$groups))) b <- fit$random_effects[[names(fit$groups)[g]]]
      if (is.null(b)) next
      b <- as.numeric(b)
      if (length(b) == ncol(Z)) random_link[ii] <- as.vector(Z[ii, , drop = FALSE] %*% b)
    }
  }
  link <- fixed_link + random_link
  family <- .memwas_normalize_family(fit$family %||% settings$family %||% "gaussian")
  mu <- if (family == "gaussian") link else {
    parts <- .memwas_family_parts(family, y = y, control = settings$control %||% list())
    parts$linkinv(link)
  }
  list(y = y, fitted_link = link, fitted_response = as.numeric(mu))
}

.memwas_extract_assumption_data <- function(fit, settings) {
  cf <- .memwas_conditional_fitted(fit, settings)
  y <- cf$y
  mu <- as.numeric(cf$fitted_response)
  eta <- as.numeric(cf$fitted_link)
  if (!length(y) || length(mu) != length(y)) stop("Fitted values or residuals are unavailable for assumption screening.", call. = FALSE)
  family <- .memwas_normalize_family(fit$family %||% settings$family %||% "gaussian")
  resid <- y - mu
  varmu <- .memwas_assumption_family_variance(family, mu, fit, settings)
  pearson <- resid / sqrt(pmax(varmu, 1e-8))
  gt <- .memwas_group_time_from_fit(fit, settings, length(y))
  list(fit = fit, settings = settings, family = family, y = y, fitted = mu,
       fitted_link = eta, residual = resid, pearson = pearson, variance = varmu,
       group = gt$group, time = gt$time, n = length(y))
}

.memwas_runs_test <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(c(statistic = NA_real_, runs = NA_real_, p.value = NA_real_))
  s <- sign(x)
  s <- s[s != 0]
  if (length(s) < 4L || length(unique(s)) < 2L) return(c(statistic = NA_real_, runs = NA_real_, p.value = NA_real_))
  n1 <- sum(s > 0)
  n2 <- sum(s < 0)
  runs <- 1 + sum(s[-1L] != s[-length(s)])
  mean_runs <- 1 + 2 * n1 * n2 / (n1 + n2)
  var_runs <- (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / (((n1 + n2)^2) * (n1 + n2 - 1))
  if (!is.finite(var_runs) || var_runs <= 0) return(c(statistic = NA_real_, runs = runs, p.value = NA_real_))
  z <- (runs - mean_runs) / sqrt(var_runs)
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  c(statistic = z, runs = runs, p.value = p)
}

.memwas_ordered_residuals <- function(e, d) {
  e <- as.numeric(e)
  if (!is.null(d$group) && !is.null(d$time) && length(d$group) == length(e) && length(d$time) == length(e)) {
    return(e[order(d$group, d$time)])
  }
  e
}

.memwas_lag_pairs <- function(e, group, time, lag = 1L) {
  e <- as.numeric(e)
  lag <- as.integer(lag[1L])
  if (!is.finite(lag) || lag < 1L) lag <- 1L
  x <- numeric(0L)
  y <- numeric(0L)
  if (!is.null(group) && length(group) == length(e)) {
    if (is.null(time) || length(time) != length(e)) time <- ave(seq_along(e), group, FUN = seq_along)
    ids <- unique(as.character(group))
    for (id in ids) {
      ii <- which(as.character(group) == id)
      ii <- ii[order(time[ii])]
      if (length(ii) > lag) {
        x <- c(x, e[ii[seq_len(length(ii) - lag)]])
        y <- c(y, e[ii[seq.int(lag + 1L, length(ii))]])
      }
    }
  } else if (length(e) > lag) {
    x <- e[seq_len(length(e) - lag)]
    y <- e[seq.int(lag + 1L, length(e))]
  }
  ok <- is.finite(x) & is.finite(y)
  list(x = x[ok], y = y[ok])
}

.memwas_lag_correlation_test <- function(e, group, time, lag = 1L) {
  pr <- .memwas_lag_pairs(e, group, time, lag = lag)
  if (length(pr$x) < 4L || stats::sd(pr$x) == 0 || stats::sd(pr$y) == 0) {
    return(c(statistic = NA_real_, correlation = NA_real_, p.value = NA_real_, n = length(pr$x)))
  }
  r <- suppressWarnings(stats::cor(pr$x, pr$y))
  if (!is.finite(r) || abs(r) >= 1) return(c(statistic = NA_real_, correlation = r, p.value = NA_real_, n = length(pr$x)))
  tstat <- r * sqrt((length(pr$x) - 2) / pmax(1 - r^2, 1e-8))
  p <- 2 * stats::pt(abs(tstat), df = length(pr$x) - 2, lower.tail = FALSE)
  c(statistic = tstat, correlation = r, p.value = p, n = length(pr$x))
}

.memwas_durbin_watson_test <- function(e, group, time) {
  pr <- .memwas_lag_pairs(e, group, time, lag = 1L)
  den <- sum(as.numeric(e)^2, na.rm = TRUE)
  dw <- if (length(pr$x) && is.finite(den) && den > 0) sum((pr$y - pr$x)^2, na.rm = TRUE) / den else NA_real_
  ct <- .memwas_lag_correlation_test(e, group, time, lag = 1L)
  c(statistic = dw, correlation = ct["correlation"], p.value = ct["p.value"], n = ct["n"])
}

.memwas_ljung_box_pooled <- function(e, group, time, max_lag = 10L) {
  e <- as.numeric(e)
  max_lag <- as.integer(max_lag[1L])
  if (!is.finite(max_lag) || max_lag < 1L) max_lag <- 10L
  if (!is.null(group) && length(group) == length(e)) {
    max_obs_lag <- max(as.integer(table(group)) - 1L, na.rm = TRUE)
  } else {
    max_obs_lag <- length(e) - 1L
  }
  h <- min(max_lag, max_obs_lag)
  if (!is.finite(h) || h < 1L) return(c(statistic = NA_real_, df = NA_real_, p.value = NA_real_, n = length(e)))
  rho <- numeric(0L)
  npair <- numeric(0L)
  for (k in seq_len(h)) {
    ct <- .memwas_lag_correlation_test(e, group, time, lag = k)
    if (is.finite(ct["correlation"])) {
      rho <- c(rho, ct["correlation"])
      npair <- c(npair, ct["n"])
    }
  }
  m <- length(rho)
  if (m == 0L) return(c(statistic = NA_real_, df = NA_real_, p.value = NA_real_, n = length(e)))
  n_eff <- max(npair, na.rm = TRUE) + 1L
  q <- n_eff * (n_eff + 2) * sum((rho^2) / pmax(n_eff - seq_len(m), 1), na.rm = TRUE)
  p <- stats::pchisq(q, df = m, lower.tail = FALSE)
  c(statistic = q, df = m, p.value = p, n = n_eff)
}

.memwas_check_autocorrelation <- function(method, d, alpha) {
  e <- d$pearson
  if (length(e) < 4L) return(.memwas_assumption_empty_result("Autocorrelation", method, "Residuals unavailable or too few observations.", alpha))
  if (method == "DurbinWatson") {
    z <- .memwas_durbin_watson_test(e, d$group, d$time)
    return(.memwas_assumption_result("Autocorrelation", method, "Durbin-Watson with pooled lag-1 p-value", "Pearson residuals", z["statistic"], 1, z["correlation"], z["p.value"], z["n"], alpha, "Auxiliary column is pooled lag-1 residual correlation; p-value uses a t approximation."))
  }
  if (method == "LjungBox") {
    max_lag <- as.integer(d$settings$control$assumption_lag %||% d$settings$control$ljung_box_lag %||% 10L)
    z <- .memwas_ljung_box_pooled(e, d$group, d$time, max_lag = max_lag)
    return(.memwas_assumption_result("Autocorrelation", method, "Pooled within-cluster Ljung-Box portmanteau", "Pearson residuals", z["statistic"], z["df"], NA_real_, z["p.value"], z["n"], alpha, "Uses within-cluster lag pairs and avoids cross-cluster residual pairs."))
  }
  if (method == "Lag1Correlation") {
    z <- .memwas_lag_correlation_test(e, d$group, d$time, lag = 1L)
    return(.memwas_assumption_result("Autocorrelation", method, "Pooled lag-1 residual correlation", "Pearson residuals", z["statistic"], z["n"] - 2, z["correlation"], z["p.value"], z["n"], alpha, "Auxiliary column is the lag-1 correlation."))
  }
  if (method == "Runs") {
    z <- .memwas_runs_test(.memwas_ordered_residuals(e, d))
    return(.memwas_assumption_result("Autocorrelation", method, "Runs test of ordered residual signs", "Pearson residuals", z["statistic"], NA_real_, z["runs"], z["p.value"], length(e), alpha, "Auxiliary column is the number of runs."))
  }
  .memwas_assumption_empty_result("Autocorrelation", method, "Unsupported autocorrelation method.", alpha)
}

.memwas_jarque_bera <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4L || stats::sd(x) == 0) return(c(statistic = NA_real_, skewness = NA_real_, excess_kurtosis = NA_real_, p.value = NA_real_, n = n))
  z <- (x - mean(x)) / stats::sd(x)
  skew <- mean(z^3)
  exk <- mean(z^4) - 3
  jb <- n / 6 * (skew^2 + exk^2 / 4)
  p <- stats::pchisq(jb, df = 2L, lower.tail = FALSE)
  c(statistic = jb, skewness = skew, excess_kurtosis = exk, p.value = p, n = n)
}

.memwas_assumption_deviance <- function(d) {
  y <- pmax(as.numeric(d$y), 0)
  mu <- pmax(as.numeric(d$fitted), 1e-8)
  fam <- d$family
  if (fam == "gaussian") {
    sig2 <- pmax(mean(d$residual^2, na.rm = TRUE), 1e-8)
    return(sum(d$residual^2 / sig2, na.rm = TRUE))
  }
  if (fam == "binomial") {
    mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
    y <- pmin(pmax(y, 0), 1)
    term1 <- ifelse(y == 0, 0, y * log(pmax(y, 1e-8) / mu))
    term2 <- ifelse(y == 1, 0, (1 - y) * log(pmax(1 - y, 1e-8) / pmax(1 - mu, 1e-8)))
    return(2 * sum(term1 + term2, na.rm = TRUE))
  }
  if (fam == "poisson") {
    term <- ifelse(y == 0, 0, y * log(pmax(y, 1e-8) / mu)) - (y - mu)
    return(2 * sum(term, na.rm = TRUE))
  }
  if (fam == "negative_binomial") {
    theta <- suppressWarnings(as.numeric(d$fit$family_parameters$theta %||% d$fit$family_parameters$size %||% d$settings$control$negative_binomial_theta %||% NA_real_)[1L])
    if (!is.finite(theta) || theta <= 0) theta <- .memwas_estimate_nb_theta(y, d$settings$control %||% list())
    term1 <- ifelse(y == 0, 0, y * log(pmax(y, 1e-8) / mu))
    term2 <- (y + theta) * log((y + theta) / pmax(mu + theta, 1e-8))
    return(2 * sum(term1 - term2, na.rm = TRUE))
  }
  if (fam %in% c("gamma", "exponential")) {
    y <- pmax(as.numeric(d$y), 1e-8)
    return(2 * sum((y - mu) / mu - log(y / mu), na.rm = TRUE))
  }
  NA_real_
}

.memwas_quantile_residuals <- function(d) {
  y <- as.numeric(d$y)
  mu <- as.numeric(d$fitted)
  fam <- d$family
  u <- rep(NA_real_, length(y))
  if (fam == "gaussian") {
    sig <- sqrt(pmax(d$variance, 1e-8))
    u <- stats::pnorm(y, mean = mu, sd = sig)
  } else if (fam == "binomial") {
    yy <- round(y)
    pp <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
    u <- stats::pbinom(yy - 1L, size = 1L, prob = pp) + 0.5 * stats::dbinom(yy, size = 1L, prob = pp)
  } else if (fam == "poisson") {
    yy <- round(pmax(y, 0))
    u <- stats::ppois(yy - 1L, lambda = pmax(mu, 1e-8)) + 0.5 * stats::dpois(yy, lambda = pmax(mu, 1e-8))
  } else if (fam == "negative_binomial") {
    yy <- round(pmax(y, 0))
    theta <- suppressWarnings(as.numeric(d$fit$family_parameters$theta %||% d$fit$family_parameters$size %||% d$settings$control$negative_binomial_theta %||% NA_real_)[1L])
    if (!is.finite(theta) || theta <= 0) theta <- .memwas_estimate_nb_theta(yy, d$settings$control %||% list())
    u <- stats::pnbinom(yy - 1L, size = theta, mu = pmax(mu, 1e-8)) + 0.5 * stats::dnbinom(yy, size = theta, mu = pmax(mu, 1e-8))
  } else if (fam == "gamma") {
    shape <- suppressWarnings(as.numeric(d$fit$family_parameters$shape %||% d$settings$control$gamma_shape %||% NA_real_)[1L])
    if (!is.finite(shape) || shape <= 0) shape <- .memwas_estimate_gamma_shape(y, d$settings$control %||% list())
    u <- stats::pgamma(pmax(y, 1e-8), shape = shape, rate = shape / pmax(mu, 1e-8))
  } else if (fam == "exponential") {
    u <- stats::pexp(pmax(y, 1e-8), rate = 1 / pmax(mu, 1e-8))
  }
  stats::qnorm(pmin(pmax(u, 1e-8), 1 - 1e-8))
}

.memwas_link_test <- function(d) {
  y <- d$y
  eta <- d$fitted_link
  ok <- is.finite(y) & is.finite(eta)
  y <- y[ok]
  eta <- eta[ok]
  if (length(y) < 8L || stats::sd(eta) == 0) return(c(statistic = NA_real_, p.value = NA_real_, n = length(y)))
  dd <- data.frame(y = y, eta = eta, eta2 = eta^2)
  fam <- d$family
  fit2 <- try(if (fam == "gaussian") {
    stats::lm(y ~ eta + eta2, data = dd)
  } else if (fam == "binomial") {
    suppressWarnings(stats::glm(y ~ eta + eta2, data = dd, family = stats::binomial()))
  } else if (fam %in% c("poisson", "negative_binomial")) {
    suppressWarnings(stats::glm(y ~ eta + eta2, data = dd, family = stats::poisson()))
  } else if (fam %in% c("gamma", "exponential")) {
    suppressWarnings(stats::glm(y ~ eta + eta2, data = dd, family = stats::Gamma(link = "log")))
  } else {
    stats::lm(y ~ eta + eta2, data = dd)
  }, silent = TRUE)
  if (inherits(fit2, "try-error")) return(c(statistic = NA_real_, p.value = NA_real_, n = length(y)))
  cc <- try(stats::coef(summary(fit2)), silent = TRUE)
  if (inherits(cc, "try-error") || is.null(dim(cc)) || !"eta2" %in% rownames(cc)) {
    return(c(statistic = NA_real_, p.value = NA_real_, n = length(y)))
  }
  pcol <- grep("Pr\\(", colnames(cc), value = TRUE)
  scol <- grep("value", colnames(cc), value = TRUE)
  p <- if (length(pcol)) cc["eta2", pcol[1L]] else NA_real_
  st <- if (length(scol)) cc["eta2", scol[1L]] else cc["eta2", ncol(cc) - 1L]
  c(statistic = st, p.value = p, n = length(y))
}

.memwas_quantile_groups <- function(x, g = 10L) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) < 6L) return(NULL)
  g <- as.integer(g[1L])
  if (!is.finite(g) || g < 2L) g <- 10L
  br <- unique(as.numeric(stats::quantile(x[ok], probs = seq(0, 1, length.out = g + 1L), na.rm = TRUE, names = FALSE, type = 7)))
  if (length(br) < 3L) return(NULL)
  out <- rep(NA_integer_, length(x))
  out[ok] <- as.integer(cut(x[ok], breaks = br, include.lowest = TRUE, labels = FALSE))
  out
}

.memwas_grouped_calibration <- function(d) {
  groups_n <- as.integer(d$settings$control$assumption_groups %||% d$settings$control$calibration_groups %||% 10L)
  g <- .memwas_quantile_groups(d$fitted, groups_n)
  if (is.null(g)) return(c(statistic = NA_real_, df = NA_real_, p.value = NA_real_, n = length(d$y)))
  lev <- sort(unique(g[!is.na(g)]))
  stat <- 0
  used <- 0L
  for (lv in lev) {
    ii <- which(g == lv)
    obs <- sum(d$y[ii], na.rm = TRUE)
    exp <- sum(d$fitted[ii], na.rm = TRUE)
    vv <- sum(d$variance[ii], na.rm = TRUE)
    if (is.finite(vv) && vv > 1e-8) {
      stat <- stat + (obs - exp)^2 / vv
      used <- used + 1L
    }
  }
  df <- max(1L, used - 2L)
  p <- if (used >= 3L) stats::pchisq(stat, df = df, lower.tail = FALSE) else NA_real_
  c(statistic = stat, df = df, p.value = p, n = length(d$y))
}

.memwas_check_distribution_link <- function(method, d, alpha) {
  n <- length(d$y)
  p <- if (!is.null(d$fit$X)) qr(d$fit$X)$rank else 1L
  df <- max(1L, n - p)
  if (method == "PearsonDispersion") {
    stat <- sum(d$pearson^2, na.rm = TRUE)
    dispersion <- stat / df
    pval <- 2 * min(stats::pchisq(stat, df = df), stats::pchisq(stat, df = df, lower.tail = FALSE))
    return(.memwas_assumption_result("Correct distribution and link function", method, "Pearson dispersion goodness-of-fit", "Pearson residuals", stat, df, dispersion, pval, n, alpha, "Auxiliary column is dispersion; two-sided chi-square p-value screens over- or under-dispersion."))
  }
  if (method == "DevianceDispersion") {
    stat <- .memwas_assumption_deviance(d)
    dispersion <- stat / df
    pval <- if (is.finite(stat)) 2 * min(stats::pchisq(stat, df = df), stats::pchisq(stat, df = df, lower.tail = FALSE)) else NA_real_
    return(.memwas_assumption_result("Correct distribution and link function", method, "Residual deviance dispersion goodness-of-fit", "Response distribution", stat, df, dispersion, pval, n, alpha, "Auxiliary column is deviance dispersion."))
  }
  if (method == "LinkTest") {
    z <- .memwas_link_test(d)
    note <- if (d$family == "negative_binomial") "Negative-binomial link screen uses a Poisson auxiliary GLM because base R has no negative-binomial GLM fitter." else "Auxiliary model adds squared fitted link term; significant eta2 suggests link or linear-predictor misspecification."
    return(.memwas_assumption_result("Correct distribution and link function", method, "Pregibon-style link test", "Squared fitted link", z["statistic"], 1, NA_real_, z["p.value"], z["n"], alpha, note))
  }
  if (method == "QuantileResidualNormality") {
    qres <- .memwas_quantile_residuals(d)
    qres <- qres[is.finite(qres)]
    if (length(qres) >= 3L && length(qres) <= 5000L) {
      st <- try(stats::shapiro.test(qres), silent = TRUE)
      if (!inherits(st, "try-error")) {
        return(.memwas_assumption_result("Correct distribution and link function", method, "Shapiro-Wilk on mid-P quantile residuals", "Quantile residuals", unname(st$statistic), NA_real_, NA_real_, st$p.value, length(qres), alpha, "Discrete families use deterministic mid-P quantile residuals to avoid simulation dependence."))
      }
    }
    jb <- .memwas_jarque_bera(qres)
    return(.memwas_assumption_result("Correct distribution and link function", method, "Jarque-Bera on mid-P quantile residuals", "Quantile residuals", jb["statistic"], 2, jb["skewness"], jb["p.value"], jb["n"], alpha, "Auxiliary column is skewness; used when Shapiro-Wilk is unavailable or sample size exceeds 5000."))
  }
  if (method == "GroupedCalibration") {
    z <- .memwas_grouped_calibration(d)
    return(.memwas_assumption_result("Correct distribution and link function", method, "Grouped observed-vs-expected calibration", "Fitted-value quantile groups", z["statistic"], z["df"], NA_real_, z["p.value"], z["n"], alpha, "Binomial models correspond to a Hosmer-Lemeshow-style screen; other families use the fitted variance function."))
  }
  .memwas_assumption_empty_result("Correct distribution and link function", method, "Unsupported distribution/link method.", alpha)
}

.memwas_check_conditional_independence <- function(method, d, alpha) {
  e <- d$pearson
  if (method == "WithinClusterLag1") {
    z <- .memwas_lag_correlation_test(e, d$group, d$time, lag = 1L)
    return(.memwas_assumption_result("Conditional independence of residuals", method, "Within-cluster lag-1 residual correlation", "Pearson residuals", z["statistic"], z["n"] - 2, z["correlation"], z["p.value"], z["n"], alpha, "Auxiliary column is lag-1 correlation after conditioning on fitted random effects."))
  }
  if (method == "ClusterMeanResidual") {
    if (is.null(d$group) || length(d$group) != length(e)) return(.memwas_assumption_empty_result("Conditional independence of residuals", method, "Group variable unavailable.", alpha))
    g <- as.factor(d$group)
    if (nlevels(g) < 2L || length(e) <= nlevels(g)) return(.memwas_assumption_empty_result("Conditional independence of residuals", method, "Too few observations per cluster for residual-by-cluster ANOVA.", alpha))
    fit <- try(stats::lm(e ~ g), silent = TRUE)
    if (inherits(fit, "try-error")) return(.memwas_assumption_empty_result("Conditional independence of residuals", method, "Residual-by-cluster ANOVA failed.", alpha))
    aa <- stats::anova(fit)
    return(.memwas_assumption_result("Conditional independence of residuals", method, "Residual mean differs by cluster", "Cluster", aa[1L, "F value"], aa[1L, "Df"], NA_real_, aa[1L, "Pr(>F)"], length(e), alpha, "Screens residual cluster structure not captured by fitted effects."))
  }
  if (method == "Runs") {
    z <- .memwas_runs_test(.memwas_ordered_residuals(e, d))
    return(.memwas_assumption_result("Conditional independence of residuals", method, "Runs test of conditional residual signs", "Pearson residuals", z["statistic"], NA_real_, z["runs"], z["p.value"], length(e), alpha, "Auxiliary column is the number of runs."))
  }
  .memwas_assumption_empty_result("Conditional independence of residuals", method, "Unsupported conditional-independence method.", alpha)
}

.memwas_random_effect_matrix <- function(fit) {
  re <- fit$random_effects %||% NULL
  if (!is.list(re) || length(re) == 0L) return(NULL)
  rows <- vector("list", length(re))
  for (i in seq_along(re)) {
    b <- as.numeric(re[[i]])
    if (!length(b)) next
    rows[[i]] <- b
  }
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(NULL)
  q <- max(vapply(rows, length, integer(1L)))
  mat <- matrix(NA_real_, nrow = length(rows), ncol = q)
  for (i in seq_along(rows)) mat[i, seq_along(rows[[i]])] <- rows[[i]]
  cn <- colnames(fit$Z %||% matrix(ncol = q))
  if (is.null(cn) || length(cn) != q) cn <- paste0("random_effect_", seq_len(q))
  colnames(mat) <- cn
  rn <- names(re)
  if (!is.null(rn) && length(rn) >= nrow(mat)) rownames(mat) <- rn[seq_len(nrow(mat))]
  mat
}

.memwas_check_re_normality <- function(method, d, alpha) {
  mat <- .memwas_random_effect_matrix(d$fit)
  if (is.null(mat) || nrow(mat) < 3L) return(.memwas_assumption_empty_result("Normality of random effects", method, "Random-effect estimates unavailable or too few clusters.", alpha))
  rows <- list()
  for (j in seq_len(ncol(mat))) {
    v <- mat[, j]
    v <- v[is.finite(v)]
    target <- colnames(mat)[j]
    if (length(v) < 3L || stats::sd(v) == 0) {
      rows[[length(rows) + 1L]] <- .memwas_assumption_empty_result("Normality of random effects", method, "Too few nonconstant empirical random effects.", alpha, target = target)
      next
    }
    if (method == "Shapiro") {
      vv <- v
      if (length(vv) > 5000L) vv <- vv[unique(round(seq(1, length(vv), length.out = 5000L)))]
      st <- try(stats::shapiro.test(vv), silent = TRUE)
      rows[[length(rows) + 1L]] <- if (inherits(st, "try-error")) {
        .memwas_assumption_empty_result("Normality of random effects", method, "Shapiro-Wilk test failed.", alpha, target = target)
      } else {
        .memwas_assumption_result("Normality of random effects", method, "Shapiro-Wilk on empirical random effects", target, unname(st$statistic), NA_real_, NA_real_, st$p.value, length(v), alpha, "BLUP/empirical random effects are shrunken estimates; interpret as a screening diagnostic.")
      }
    } else if (method == "JarqueBera") {
      jb <- .memwas_jarque_bera(v)
      rows[[length(rows) + 1L]] <- .memwas_assumption_result("Normality of random effects", method, "Jarque-Bera on empirical random effects", target, jb["statistic"], 2, jb["skewness"], jb["p.value"], jb["n"], alpha, "Auxiliary column is skewness; BLUP/empirical random effects are shrunken estimates.")
    } else if (method == "SkewKurtosis") {
      jb <- .memwas_jarque_bera(v)
      rows[[length(rows) + 1L]] <- .memwas_assumption_result("Normality of random effects", method, "Skewness/excess-kurtosis screen", target, jb["skewness"], NA_real_, jb["excess_kurtosis"], jb["p.value"], jb["n"], alpha, "Statistic is skewness; auxiliary column is excess kurtosis; p-value uses Jarque-Bera approximation.")
    } else {
      rows[[length(rows) + 1L]] <- .memwas_assumption_empty_result("Normality of random effects", method, "Unsupported random-effects-normality method.", alpha, target = target)
    }
  }
  do.call(rbind, rows)
}

.memwas_group_mean_matrix <- function(fit) {
  X <- fit$X %||% NULL
  groups <- fit$groups %||% NULL
  if (is.null(X) || is.null(groups) || ncol(X) < 2L || length(groups) == 0L) return(NULL)
  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]
  if (ncol(X) == 0L) return(NULL)
  keep2 <- apply(X, 2L, function(z) stats::sd(z, na.rm = TRUE) > 1e-10)
  X <- X[, keep2, drop = FALSE]
  if (ncol(X) == 0L) return(NULL)
  out <- matrix(NA_real_, nrow = length(groups), ncol = ncol(X))
  rownames(out) <- names(groups) %||% as.character(seq_along(groups))
  colnames(out) <- colnames(X)
  for (g in seq_along(groups)) out[g, ] <- colMeans(X[groups[[g]], , drop = FALSE], na.rm = TRUE)
  out
}

.memwas_correlation_min_test <- function(y, X, method = "pearson") {
  pvals <- rep(NA_real_, ncol(X))
  stats_val <- rep(NA_real_, ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    ok <- is.finite(y) & is.finite(x)
    if (sum(ok) < 4L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) next
    tt <- try(suppressWarnings(stats::cor.test(y[ok], x[ok], method = method, exact = FALSE)), silent = TRUE)
    if (!inherits(tt, "try-error")) {
      pvals[j] <- tt$p.value
      stats_val[j] <- as.numeric(tt$estimate)
    }
  }
  if (all(!is.finite(pvals))) return(list(statistic = NA_real_, p.value = NA_real_, predictor = NA_character_))
  j <- which.min(pvals)
  list(statistic = stats_val[j], p.value = pvals[j], predictor = colnames(X)[j])
}

.memwas_check_re_predictor_independence <- function(method, d, alpha) {
  re <- .memwas_random_effect_matrix(d$fit)
  gm <- .memwas_group_mean_matrix(d$fit)
  if (is.null(re) || is.null(gm)) return(.memwas_assumption_empty_result("Independence between random effects and predictors", method, "Random effects or fixed-effect predictor group means unavailable.", alpha))
  common <- intersect(rownames(re), rownames(gm))
  if (length(common) < 5L) return(.memwas_assumption_empty_result("Independence between random effects and predictors", method, "Too few clusters matched between random effects and predictor group means.", alpha))
  re <- re[common, , drop = FALSE]
  gm <- gm[common, , drop = FALSE]
  rows <- list()
  for (j in seq_len(ncol(re))) {
    y <- re[, j]
    target <- colnames(re)[j]
    if (method == "GroupMeanAssociation") {
      dd <- data.frame(.re = y, gm, check.names = FALSE)
      fit <- try(stats::lm(.re ~ ., data = dd), silent = TRUE)
      if (inherits(fit, "try-error") || is.null(summary(fit)$fstatistic)) {
        rows[[length(rows) + 1L]] <- .memwas_assumption_empty_result("Independence between random effects and predictors", method, "Mundlak-style auxiliary regression failed or was rank deficient.", alpha, target = target)
      } else {
        fs <- summary(fit)$fstatistic
        p <- stats::pf(fs[1L], fs[2L], fs[3L], lower.tail = FALSE)
        rows[[length(rows) + 1L]] <- .memwas_assumption_result("Independence between random effects and predictors", method, "Mundlak-style RE ~ group predictor means", target, fs[1L], fs[2L], ncol(gm), p, length(y), alpha, "Auxiliary column is the number of predictor group-mean columns in the screen.")
      }
    } else if (method == "Correlation") {
      ct <- .memwas_correlation_min_test(y, gm, method = "pearson")
      rows[[length(rows) + 1L]] <- .memwas_assumption_result("Independence between random effects and predictors", method, "Minimum Pearson RE-predictor mean correlation p-value", target, ct$statistic, NA_real_, NA_real_, ct$p.value, length(y), alpha, paste0("Screened predictor mean: ", ct$predictor))
    } else if (method == "RankCorrelation") {
      ct <- .memwas_correlation_min_test(y, gm, method = "spearman")
      rows[[length(rows) + 1L]] <- .memwas_assumption_result("Independence between random effects and predictors", method, "Minimum Spearman RE-predictor mean correlation p-value", target, ct$statistic, NA_real_, NA_real_, ct$p.value, length(y), alpha, paste0("Screened predictor mean: ", ct$predictor))
    } else {
      rows[[length(rows) + 1L]] <- .memwas_assumption_empty_result("Independence between random effects and predictors", method, "Unsupported RE-predictor-independence method.", alpha, target = target)
    }
  }
  do.call(rbind, rows)
}

.memwas_levene_test <- function(y, g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]
  g <- as.factor(g[ok])
  tab <- table(g)
  if (length(y) < 6L || sum(tab >= 2L) < 2L) return(c(statistic = NA_real_, df = NA_real_, p.value = NA_real_))
  keep <- names(tab)[tab >= 1L]
  sel <- g %in% keep
  y <- y[sel]
  g <- droplevels(g[sel])
  med <- tapply(y, g, stats::median, na.rm = TRUE)
  z <- abs(y - med[as.character(g)])
  fit <- try(stats::lm(z ~ g), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(statistic = NA_real_, df = NA_real_, p.value = NA_real_))
  aa <- stats::anova(fit)
  c(statistic = aa[1L, "F value"], df = aa[1L, "Df"], p.value = aa[1L, "Pr(>F)"])
}

.memwas_check_variance <- function(method, d, alpha) {
  e <- as.numeric(d$pearson)
  fv <- as.numeric(d$fitted)
  ok <- is.finite(e) & is.finite(fv)
  if (method == "BreuschPagan") {
    if (sum(ok) < 6L || stats::sd(fv[ok]) == 0) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "Fitted values unavailable or invariant.", alpha))
    aux <- try(stats::lm(I(e[ok]^2) ~ fv[ok]), silent = TRUE)
    if (inherits(aux, "try-error")) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "Breusch-Pagan auxiliary regression failed.", alpha))
    r2 <- summary(aux)$r.squared
    stat <- sum(ok) * r2
    p <- stats::pchisq(stat, df = 1L, lower.tail = FALSE)
    return(.memwas_assumption_result("Homogeneity of variance", method, "Breusch-Pagan fitted-value score", "Pearson residuals", stat, 1, r2, p, sum(ok), alpha, "Auxiliary column is R-squared from squared-residual regression."))
  }
  if (method == "White") {
    if (sum(ok) < 8L || stats::sd(fv[ok]) == 0) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "Fitted values unavailable or invariant.", alpha))
    aux <- try(stats::lm(I(e[ok]^2) ~ fv[ok] + I(fv[ok]^2)), silent = TRUE)
    if (inherits(aux, "try-error")) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "White auxiliary regression failed.", alpha))
    r2 <- summary(aux)$r.squared
    stat <- sum(ok) * r2
    p <- stats::pchisq(stat, df = 2L, lower.tail = FALSE)
    return(.memwas_assumption_result("Homogeneity of variance", method, "White fitted-value quadratic test", "Pearson residuals", stat, 2, r2, p, sum(ok), alpha, "Auxiliary column is R-squared from squared-residual regression."))
  }
  if (method == "LeveneFitted") {
    if (sum(ok) < 8L) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "Fitted values unavailable.", alpha))
    g <- .memwas_quantile_groups(fv[ok], g = 4L)
    if (is.null(g)) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "Could not form fitted-value quartile groups.", alpha))
    z <- .memwas_levene_test(e[ok], g)
    return(.memwas_assumption_result("Homogeneity of variance", method, "Brown-Forsythe/Levene test across fitted quartiles", "Fitted-value quartiles", z["statistic"], z["df"], NA_real_, z["p.value"], sum(ok), alpha, "Uses absolute deviations from group medians."))
  }
  if (method == "LeveneGroup") {
    if (is.null(d$group) || length(d$group) != length(e)) return(.memwas_assumption_empty_result("Homogeneity of variance", method, "Group variable unavailable.", alpha))
    z <- .memwas_levene_test(e, d$group)
    return(.memwas_assumption_result("Homogeneity of variance", method, "Brown-Forsythe/Levene test across clusters", "Cluster", z["statistic"], z["df"], NA_real_, z["p.value"], length(e), alpha, "Uses absolute deviations from cluster medians."))
  }
  .memwas_assumption_empty_result("Homogeneity of variance", method, "Unsupported homogeneity-of-variance method.", alpha)
}

.memwas_assumption_unavailable_table <- function(spec, reason, alpha = 0.05) {
  if (is.null(spec) || .memwas_assumption_is_empty_spec(spec)) return(.memwas_assumption_empty_table())
  labels <- c(autocorrelation_check = "Autocorrelation",
              distribution_link_check = "Correct distribution and link function",
              conditional_independence_check = "Conditional independence of residuals",
              random_effects_normality_check = "Normality of random effects",
              random_effects_predictor_independence_check = "Independence between random effects and predictors",
              homogeneity_variance_check = "Homogeneity of variance")
  rows <- list()
  for (nm in names(labels)) {
    for (method in spec[[nm]] %||% character(0L)) {
      rows[[length(rows) + 1L]] <- .memwas_assumption_empty_result(labels[[nm]], method, reason, alpha)
    }
  }
  if (length(rows)) do.call(rbind, rows) else .memwas_assumption_empty_table()
}

.memwas_run_assumption_checks <- function(fit, settings, spec = NULL, stage = "setup") {
  alpha <- .memwas_assumption_alpha(settings)
  if (is.null(spec)) spec <- settings$assumption_check_spec %||% NULL
  if (is.null(spec) || .memwas_assumption_is_empty_spec(spec)) {
    tab <- .memwas_assumption_empty_table()
    return(list(stage = stage, alpha = alpha, methods = spec, summary = tab,
                tables = .memwas_split_assumption_table(tab), note = "No GLMM assumption screening methods were requested."))
  }
  d <- try(.memwas_extract_assumption_data(fit, settings), silent = TRUE)
  if (inherits(d, "try-error")) {
    tab <- .memwas_assumption_unavailable_table(spec, as.character(d), alpha)
    return(list(stage = stage, alpha = alpha, methods = spec, summary = tab,
                tables = .memwas_split_assumption_table(tab), note = "Assumption screening could not extract residual information from the fitted model."))
  }
  rows <- list()
  for (method in spec$autocorrelation_check %||% character(0L)) {
    rows[[length(rows) + 1L]] <- .memwas_check_autocorrelation(method, d, alpha)
  }
  for (method in spec$distribution_link_check %||% character(0L)) {
    rows[[length(rows) + 1L]] <- .memwas_check_distribution_link(method, d, alpha)
  }
  for (method in spec$conditional_independence_check %||% character(0L)) {
    rows[[length(rows) + 1L]] <- .memwas_check_conditional_independence(method, d, alpha)
  }
  for (method in spec$random_effects_normality_check %||% character(0L)) {
    rows[[length(rows) + 1L]] <- .memwas_check_re_normality(method, d, alpha)
  }
  for (method in spec$random_effects_predictor_independence_check %||% character(0L)) {
    rows[[length(rows) + 1L]] <- .memwas_check_re_predictor_independence(method, d, alpha)
  }
  for (method in spec$homogeneity_variance_check %||% character(0L)) {
    rows[[length(rows) + 1L]] <- .memwas_check_variance(method, d, alpha)
  }
  tab <- if (length(rows)) do.call(rbind, rows) else .memwas_assumption_empty_table()
  row.names(tab) <- NULL
  list(stage = stage, alpha = alpha, methods = spec, summary = tab,
       tables = .memwas_split_assumption_table(tab),
       note = paste("Base-R GLMM assumption screening was computed from conditional fitted values and Pearson residuals at the", stage, "stage."))
}


# Backward-compatible internal aliases used by set_MEMWAS() and summary() wrappers.
.memwas_validate_assumption_checks <- .memwas_validate_assumption_check_settings
.memwas_empty_assumption_table <- .memwas_assumption_empty_table

.memwas_no_assumption_check_settings <- function() {
  .memwas_validate_assumption_check_settings(
    autocorrelation_check = "None",
    distribution_link_check = "None",
    conditional_independence_check = "None",
    random_effects_normality_check = "None",
    random_effects_predictor_independence_check = "None",
    homogeneity_variance_check = "None"
  )
}

.memwas_assumption_spec_table <- function(spec) {
  if (is.null(spec) || !is.list(spec)) spec <- .memwas_no_assumption_check_settings()
  labels <- c(
    autocorrelation_check = "Autocorrelation",
    distribution_link_check = "Correct distribution and link function",
    conditional_independence_check = "Conditional independence of residuals",
    random_effects_normality_check = "Normality of random effects",
    random_effects_predictor_independence_check = "Independence between random effects and predictors",
    homogeneity_variance_check = "Homogeneity of variance"
  )
  rows <- vector("list", length(labels))
  i <- 0L
  for (nm in names(labels)) {
    i <- i + 1L
    methods <- spec[[nm]] %||% character(0L)
    rows[[i]] <- data.frame(
      assumption = unname(labels[[nm]]),
      argument = nm,
      methods = if (length(methods)) paste(methods, collapse = "; ") else "None",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  do.call(rbind, rows)
}

.memwas_get_assumption_summary_table <- function(checks) {
  if (is.null(checks)) return(.memwas_assumption_empty_table())
  if (is.data.frame(checks)) return(checks)
  if (is.list(checks) && is.data.frame(checks$summary)) return(checks$summary)
  .memwas_assumption_empty_table()
}

.memwas_get_assumption_split_tables <- function(checks) {
  if (is.null(checks)) return(list())
  if (is.list(checks) && is.list(checks$tables)) return(checks$tables)
  .memwas_split_assumption_table(.memwas_get_assumption_summary_table(checks))
}

.memwas_assumption_checks_unavailable <- function(spec, settings, reason, stage = "setup") {
  alpha <- .memwas_assumption_alpha(settings %||% list())
  tab <- .memwas_assumption_unavailable_table(spec, reason, alpha = alpha)
  list(
    stage = stage,
    alpha = alpha,
    methods = spec,
    summary = tab,
    tables = .memwas_split_assumption_table(tab),
    note = reason
  )
}

.memwas_print_assumption_results <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L) {
    cat("No GLMM assumption-screening results available.\n")
    return(invisible(NULL))
  }
  z <- tab
  num_cols <- intersect(c("statistic", "df", "auxiliary", "p_value", "alpha"), names(z))
  for (cc in num_cols) z[[cc]] <- signif(z[[cc]], 5)
  print(z, row.names = FALSE)
  invisible(z)
}

.memwas_assumption_spec_from_settings <- function(settings, default = "None") {
  if (!is.null(settings) && !is.null(settings$assumption_check_spec)) return(settings$assumption_check_spec)
  settings <- settings %||% list()
  .memwas_validate_assumption_check_settings(
    autocorrelation_check = settings$autocorrelation_check %||% default,
    distribution_link_check = settings$distribution_link_check %||% default,
    conditional_independence_check = settings$conditional_independence_check %||% default,
    random_effects_normality_check = settings$random_effects_normality_check %||% default,
    random_effects_predictor_independence_check = settings$random_effects_predictor_independence_check %||% default,
    homogeneity_variance_check = settings$homogeneity_variance_check %||% default
  )
}
