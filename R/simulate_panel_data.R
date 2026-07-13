#' Simulate MEMWAS panel data
#'
#' @description
#' Vectorized panel-data simulation for MEMWAS-style experiments with correlated
#' predictors, predictor-specific temporal autocorrelation structures, optional
#' nonlinear fixed-effect components, subject-specific random intercepts,
#' residual autocorrelation, and Gaussian, binomial, or Poisson outcomes.
#'
#' @param n_id Integer scalar. Number of individuals.
#' @param n_time Integer scalar. Number of repeated time points per individual.
#' @param beta Numeric vector. Linear fixed-effect coefficients for predictors.
#'   The length of `beta` defines the number of predictors. If `beta` is named,
#'   the names are used as predictor names in the returned data set and in
#'   `predictor_autocor` matching.
#' @param cor_matrix Numeric matrix. Target same-time predictor correlation
#'   matrix with dimension `length(beta)` by `length(beta)`. If both `beta` and
#'   `cor_matrix` are named, `cor_matrix` is reordered to match `names(beta)`.
#' @param intercept Numeric scalar. Fixed intercept on the linear predictor scale.
#' @param sigma_eps Numeric scalar. Standard deviation of the latent residual
#'   process. For Gaussian outcomes this is the additive residual standard
#'   deviation. For binomial and Poisson outcomes this is a latent-scale
#'   Gaussian residual standard deviation.
#' @param sigma_b Numeric scalar. Standard deviation of subject-specific random
#'   intercepts.
#' @param family Character. Outcome family: `"gaussian"`, `"binomial"`, or
#'   `"poisson"`.
#' @param autocor Character. Backward-compatible residual autocorrelation
#'   structure. Supported values are `"NONE"`, `"AR1"`, `"ARp"`,
#'   `"ARMA11"`, `"CS"`, `"TOEP"`, and `"UN"`. The aliases `"AR(p)"`,
#'   `"ARMA(1,1)"`, `"TOEPLITZ"`, `"EXCHANGEABLE"`, and `"USER"` are also
#'   accepted. Ignored when `residual_autocor` is supplied.
#' @param autocor_param List. Backward-compatible residual autocorrelation
#'   parameters. Use `list(rho = ...)` for `"AR1"` and `"CS"`,
#'   `list(phi = c(...))` or `list(ar = c(...))` for `"ARp"`,
#'   `list(ar = ..., ma = ...)` or `list(phi = ..., theta = ...)` for
#'   `"ARMA11"`, `list(rho_vec = c(...))` for `"TOEP"`, and
#'   `list(R = matrix(...))` or `list(Sigma = matrix(...))` for `"UN"`.
#' @param nonlinear_fun Function or `NULL`. Optional user-defined nonlinear
#'   component. The function must take the predictor matrix and return a numeric
#'   vector of length `n_id * n_time`. If supplied, it overrides
#'   `nonlinear_type`.
#' @param nonlinear_type Character. Built-in nonlinear form: `"none"`,
#'   `"quadratic"`, or `"interaction"`.
#' @param beta_nl Numeric vector or `NULL`. Coefficients for built-in nonlinear
#'   terms. For `nonlinear_type = "quadratic"`, the length must equal
#'   `length(beta)`. For `nonlinear_type = "interaction"`, the length must equal
#'   `choose(length(beta), 2)`.
#' @param binomial_size Integer scalar. Number of trials for binomial outcomes.
#'   The default value of 1 gives Bernoulli outcomes.
#' @param exact_predictor_cor Logical scalar. If `TRUE`, the simulator attempts
#'   to make the finite-sample predictor correlation match `cor_matrix`. With no
#'   predictor autocorrelation this is applied to the stacked long-format
#'   predictor matrix, matching the previous implementation. With
#'   predictor-specific autocorrelation this is applied to the subject-level wide
#'   predictor covariance and requires `n_id > length(beta) * n_time`; otherwise
#'   population-target simulation is used with a warning.
#' @param return_latent Logical scalar. If `TRUE`, include `eta` and
#'   `latent_eps` in the returned data frame.
#' @param pd_tol Numeric scalar. Positive-definiteness tolerance.
#' @param predictor_autocor Character, list, or `NULL`. Predictor-level temporal
#'   autocorrelation specification. `NULL` means no predictor autocorrelation.
#'   A global specification can be written as `list(structure = "AR1", rho = ...)`.
#'   A predictor-specific specification can be written as
#'   `list(X1 = list(structure = "AR1", rho = 0.7),
#'         X2 = list(structure = "CS", rho = 0.3),
#'         default = list(structure = "NONE"))`.
#'   Supported structures are `"NONE"`, `"AR1"`, `"ARp"`, `"ARMA11"`,
#'   `"CS"`, `"TOEP"`, `"UN"`, and `"USER"`.
#' @param residual_autocor Character, list, or `NULL`. New explicit residual
#'   autocorrelation specification. If supplied, it overrides `autocor` and
#'   `autocor_param`. Uses the same structures and parameter names as
#'   `predictor_autocor`, but represents residual/latent-error autocorrelation.
#' @param predictor_means Numeric vector or `NULL`. Optional predictor means.
#'   If `NULL`, all predictor means are zero. Named vectors are reordered to
#'   match predictor names.
#' @param predictor_sds Numeric vector or `NULL`. Optional predictor standard
#'   deviations. If `NULL`, all predictor standard deviations are one. Named
#'   vectors are reordered to match predictor names.
#' @param cross_lag_rule Character. Rule for cross-predictor, cross-time
#'   covariance blocks when predictors have different autocorrelation structures.
#'   `"same_time_only"` imposes `cor_matrix` only at the same time point and is
#'   the safest default. `"geometric"` and `"average"` extend cross-predictor
#'   covariance across lags using the predictor autocorrelation matrices.
#'   `"zero"` sets all cross-predictor covariance blocks to zero.
#' @param make_pd Logical scalar. If `TRUE`, near or non-positive-definite
#'   predictor autocorrelation/covariance matrices are repaired by eigenvalue
#'   clipping. If `FALSE`, invalid matrices trigger an error.
#' @param return_components Logical scalar. If `TRUE`, return a list containing
#'   the simulated data and the matrices/specifications used to generate it.
#' @param phi Numeric vector or `NULL`. Optional backward-compatible shortcut for
#'   predictor-level AR(1) coefficients. If `predictor_autocor` is `NULL` and
#'   `phi` is supplied, it is converted to
#'   `predictor_autocor = list(structure = "AR1", rho = phi)`.
#'
#' @returns
#' If `return_components = FALSE`, a data frame containing `id`, `time`, `y`, and
#' predictors. If `return_latent = TRUE`, the data frame also contains `eta` and
#' `latent_eps`.
#'
#' If `return_components = TRUE`, a list containing the data frame, the
#' predictor-specific autocorrelation specifications, predictor autocorrelation
#' matrices, residual autocorrelation specification, residual autocorrelation
#' matrix, joint predictor covariance matrix, latent residuals, linear predictor,
#' and random intercepts.
#'
#' @details
#' The outcome is generated on a latent linear predictor scale:
#'
#' \deqn{\eta_{it} = \alpha + X_{it}^{\top}\beta + g(X_{it}) + b_i + e_{it},}
#'
#' where `intercept` is \eqn{\alpha}, \eqn{g(X_{it})} is the optional nonlinear
#' term, \eqn{b_i \sim N(0, \sigma_b^2)} is a subject-specific random intercept,
#' and \eqn{e_i = (e_{i1}, \ldots, e_{iT})} has covariance
#' \eqn{\sigma_\epsilon^2 R_\epsilon}.
#'
#' Predictor-specific autocorrelation is modeled separately from residual
#' autocorrelation. For predictor \eqn{j}, the within-subject predictor
#' correlation can be written as
#'
#' \deqn{\mathrm{Corr}(X_{ij t}, X_{ij s}) = R_j(t, s).}
#'
#' Residual autocorrelation is represented separately as
#'
#' \deqn{\mathrm{Corr}(e_{it}, e_{is}) = R_\epsilon(t, s).}
#'
#' This distinction lets one predictor follow an AR(1) structure, another follow
#' compound symmetry, another follow a Toeplitz or user-specified structure, while
#' the outcome residual process uses its own structure.
#'
#' The inherited arguments `autocor` and `autocor_param` continue to refer to the
#' residual/latent-error autocorrelation. The newer `residual_autocor` argument is
#' preferred when both predictor and residual autocorrelation are being specified
#' because it makes the distinction explicit.
#'
#' @examples
#' \dontrun{
#' beta <- c(X1 = 0.5, X2 = -0.3, X3 = 0.2)
#' cor_matrix <- matrix(
#'   c(1.0, 0.3, 0.1,
#'     0.3, 1.0, 0.2,
#'     0.1, 0.2, 1.0),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(names(beta), names(beta))
#' )
#'
#' dat <- simulate_panel_data(
#'   n_id = 500,
#'   n_time = 5,
#'   beta = beta,
#'   cor_matrix = cor_matrix,
#'   family = "gaussian",
#'   predictor_autocor = list(
#'     X1 = list(structure = "AR1", rho = 0.7),
#'     X2 = list(structure = "CS", rho = 0.3),
#'     X3 = list(structure = "TOEP", rho_vec = c(0.5, 0.3, 0.1, 0.0))
#'   ),
#'   residual_autocor = list(structure = "AR1", rho = 0.4),
#'   return_latent = TRUE
#' )
#'
#' head(dat)
#' }
#'
#' @noRd
.simulate_panel_data <- function(
    n_id = 1000,
    n_time = 3,
    beta,
    cor_matrix,
    intercept = 0,
    sigma_eps = 1,
    sigma_b = 1,
    family = c("gaussian", "binomial", "poisson"),
    autocor = "NONE",
    autocor_param = list(),
    nonlinear_fun = NULL,
    nonlinear_type = c("none", "quadratic", "interaction"),
    beta_nl = NULL,
    binomial_size = 1L,
    exact_predictor_cor = TRUE,
    return_latent = FALSE,
    pd_tol = 1e-8,
    predictor_autocor = NULL,
    residual_autocor = NULL,
    predictor_means = NULL,
    predictor_sds = NULL,
    cross_lag_rule = c("same_time_only", "geometric", "average", "zero"),
    make_pd = TRUE,
    return_components = FALSE,
    phi = NULL
) {

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 1. Input validation -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  family <- match.arg(family)
  nonlinear_type <- match.arg(nonlinear_type)
  cross_lag_rule <- match.arg(cross_lag_rule)

  `%||%` <- function(x, y) if (is.null(x)) y else x

  if (missing(beta)) {
    stop("beta must be supplied.", call. = FALSE)
  }
  if (missing(cor_matrix)) {
    stop("cor_matrix must be supplied.", call. = FALSE)
  }

  if (!is.numeric(n_id) || length(n_id) != 1L || !is.finite(n_id) || n_id < 1L) {
    stop("n_id must be a positive integer-like scalar.", call. = FALSE)
  }
  if (!is.numeric(n_time) || length(n_time) != 1L || !is.finite(n_time) || n_time < 1L) {
    stop("n_time must be a positive integer-like scalar.", call. = FALSE)
  }
  n_id <- as.integer(n_id)
  n_time <- as.integer(n_time)

  if (!is.numeric(beta) || length(beta) < 1L || any(!is.finite(beta))) {
    stop("beta must be a non-empty finite numeric vector.", call. = FALSE)
  }
  p <- length(beta)
  n_obs <- n_id * n_time

  cor_matrix <- as.matrix(cor_matrix)
  if (!all(dim(cor_matrix) == c(p, p))) {
    stop("cor_matrix must have dimension length(beta) by length(beta).", call. = FALSE)
  }

  beta_names <- names(beta)
  beta_named <- !is.null(beta_names) && length(beta_names) == p && all(nzchar(beta_names))
  cor_names <- colnames(cor_matrix)
  cor_named <- !is.null(cor_names) && length(cor_names) == p && all(nzchar(cor_names))

  if (beta_named) {
    var_names <- beta_names
    if (cor_named) {
      if (!all(var_names %in% cor_names) || is.null(rownames(cor_matrix)) || !all(var_names %in% rownames(cor_matrix))) {
        stop("When beta and cor_matrix are named, cor_matrix must include all beta names as row and column names.", call. = FALSE)
      }
      cor_matrix <- cor_matrix[var_names, var_names, drop = FALSE]
    }
  } else if (cor_named) {
    var_names <- cor_names
  } else {
    var_names <- paste0("X", seq_len(p))
  }

  if (!is.numeric(intercept) || length(intercept) != 1L || !is.finite(intercept)) {
    stop("intercept must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(sigma_eps) || length(sigma_eps) != 1L || !is.finite(sigma_eps) || sigma_eps < 0) {
    stop("sigma_eps must be a non-negative finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(sigma_b) || length(sigma_b) != 1L || !is.finite(sigma_b) || sigma_b < 0) {
    stop("sigma_b must be a non-negative finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(binomial_size) || length(binomial_size) != 1L || !is.finite(binomial_size) || binomial_size < 1L) {
    stop("binomial_size must be a positive integer-like scalar.", call. = FALSE)
  }
  binomial_size <- as.integer(binomial_size)

  if (!is.logical(exact_predictor_cor) || length(exact_predictor_cor) != 1L || is.na(exact_predictor_cor)) {
    stop("exact_predictor_cor must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_latent) || length(return_latent) != 1L || is.na(return_latent)) {
    stop("return_latent must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_components) || length(return_components) != 1L || is.na(return_components)) {
    stop("return_components must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(make_pd) || length(make_pd) != 1L || is.na(make_pd)) {
    stop("make_pd must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(pd_tol) || length(pd_tol) != 1L || !is.finite(pd_tol) || pd_tol <= 0) {
    stop("pd_tol must be a positive finite numeric scalar.", call. = FALSE)
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 2. Set helper functions -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  normalize_numeric_vector <- function(x, default, name, positive = FALSE) {
    if (is.null(x)) {
      out <- rep(default, p)
      names(out) <- var_names
      return(out)
    }
    if (!is.numeric(x)) {
      stop(name, " must be numeric.", call. = FALSE)
    }
    if (!is.null(names(x)) && all(var_names %in% names(x))) {
      out <- as.numeric(x[var_names])
    } else if (length(x) == 1L) {
      out <- rep(as.numeric(x), p)
    } else if (length(x) == p) {
      out <- as.numeric(x)
    } else {
      stop(name, " must have length 1, length(beta), or names matching the predictors.", call. = FALSE)
    }
    if (any(!is.finite(out))) {
      stop(name, " must contain finite values.", call. = FALSE)
    }
    if (positive && any(out <= 0)) {
      stop(name, " must contain positive values.", call. = FALSE)
    }
    names(out) <- var_names
    out
  }

  predictor_means <- normalize_numeric_vector(predictor_means, default = 0, name = "predictor_means")
  predictor_sds <- normalize_numeric_vector(predictor_sds, default = 1, name = "predictor_sds", positive = TRUE)

  normalize_autocor <- function(x) {
    if (is.null(x)) x <- "NONE"
    x <- toupper(as.character(x[1L]))
    x <- gsub("\\s+", "", x)
    x <- gsub("[-_]+", "", x)

    if (x %in% c("NONE", "NO", "INDEPENDENT", "INDEPENDENCE", "IDENTITY", "I")) return("NONE")
    if (x %in% c("AR1", "ARONE")) return("AR1")
    if (x %in% c("ARP", "AR(P)", "AR")) return("ARP")
    if (x %in% c("ARMA11", "ARMA(1,1)", "ARMA")) return("ARMA11")
    if (x %in% c("CS", "EXCHANGEABLE", "COMPOUNDSYMMETRY")) return("CS")
    if (x %in% c("TOEP", "TOEPLITZ")) return("TOEP")
    if (x %in% c("UN", "UNSTRUCTURED")) return("UN")
    if (x %in% c("USER", "CUSTOM")) return("USER")

    stop("Unsupported autocorrelation structure: ", x,
         ". Use 'NONE', 'AR1', 'ARp', 'ARMA11', 'CS', 'TOEP', 'UN', or 'USER'.",
         call. = FALSE)
  }

  check_correlation_matrix <- function(R, name, expected_dim = NULL) {
    if (!is.matrix(R) || !is.numeric(R)) {
      stop(name, " must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(R) != ncol(R)) {
      stop(name, " must be square.", call. = FALSE)
    }
    if (!is.null(expected_dim) && nrow(R) != expected_dim) {
      stop(name, " has incorrect dimension.", call. = FALSE)
    }
    if (any(!is.finite(R))) {
      stop(name, " must contain only finite values.", call. = FALSE)
    }
    if (max(abs(R - t(R))) > sqrt(.Machine$double.eps)) {
      stop(name, " must be symmetric.", call. = FALSE)
    }
    if (max(abs(diag(R) - 1)) > 1e-7) {
      stop(name, " must have diagonal elements equal to 1.", call. = FALSE)
    }
    if (any(abs(R) > 1 + 1e-10)) {
      stop(name, " must contain correlations in [-1, 1].", call. = FALSE)
    }
    ev <- eigen((R + t(R)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= pd_tol) {
      stop(name, " must be positive definite. Minimum eigenvalue = ",
           signif(min(ev), 4), ".", call. = FALSE)
    }
    invisible(TRUE)
  }

  cor_matrix <- (cor_matrix + t(cor_matrix)) / 2
  diag(cor_matrix) <- 1
  dimnames(cor_matrix) <- list(var_names, var_names)
  check_correlation_matrix(cor_matrix, "cor_matrix", expected_dim = p)

  as_positive_definite_corr <- function(R, name) {
    R <- as.matrix(R)
    if (!all(dim(R) == c(n_time, n_time))) {
      stop(name, " must have dimension n_time by n_time.", call. = FALSE)
    }
    if (any(!is.finite(R))) {
      stop(name, " must contain only finite values.", call. = FALSE)
    }
    R <- (R + t(R)) / 2
    diag(R) <- 1
    if (any(abs(R) > 1 + 1e-8)) {
      stop(name, " must contain correlations in [-1, 1].", call. = FALSE)
    }

    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) > pd_tol) {
      return(R)
    }

    if (!make_pd) {
      stop(name, " is not positive definite. Minimum eigenvalue = ",
           signif(min(ev), 4), ".", call. = FALSE)
    }

    warning(name, " was not positive definite and has been repaired by eigenvalue clipping.",
            call. = FALSE)
    eig <- eigen(R, symmetric = TRUE)
    vals <- pmax(eig$values, pd_tol)
    R2 <- eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
    R2 <- (R2 + t(R2)) / 2
    d <- sqrt(pmax(diag(R2), .Machine$double.eps))
    R2 <- sweep(sweep(R2, 1L, d, "/"), 2L, d, "/")
    diag(R2) <- 1
    R2
  }

  as_positive_definite_cov <- function(S, name) {
    S <- as.matrix(S)
    if (nrow(S) != ncol(S)) {
      stop(name, " must be square.", call. = FALSE)
    }
    if (any(!is.finite(S))) {
      stop(name, " must contain only finite values.", call. = FALSE)
    }
    S <- (S + t(S)) / 2
    target_diag <- diag(S)
    if (any(target_diag <= 0)) {
      stop(name, " must have positive diagonal values.", call. = FALSE)
    }
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) > pd_tol) {
      return(S)
    }

    if (!make_pd) {
      stop(name, " is not positive definite. Minimum eigenvalue = ",
           signif(min(ev), 4), ".", call. = FALSE)
    }

    warning(name, " was not positive definite and has been repaired by eigenvalue clipping.",
            call. = FALSE)
    eig <- eigen(S, symmetric = TRUE)
    vals <- pmax(eig$values, pd_tol)
    S2 <- eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
    S2 <- (S2 + t(S2)) / 2
    scale <- sqrt(pmax(target_diag, .Machine$double.eps) / pmax(diag(S2), .Machine$double.eps))
    S2 <- sweep(sweep(S2, 1L, scale, "*"), 2L, scale, "*")
    diag(S2) <- target_diag
    S2
  }

  get_param <- function(x, names, required = TRUE, default = NULL, label = "autocorrelation specification") {
    for (nm in names) {
      if (!is.null(x[[nm]])) return(x[[nm]])
    }
    if (required) {
      stop(label, " is missing a parameter. Expected one of: ",
           paste(names, collapse = ", "), call. = FALSE)
    }
    default
  }

  check_ar_stationarity <- function(ar, label) {
    if (!is.numeric(ar) || length(ar) < 1L || any(!is.finite(ar))) {
      stop(label, ": AR coefficients must be finite numeric values.", call. = FALSE)
    }
    roots <- polyroot(c(1, -as.numeric(ar)))
    if (any(Mod(roots) <= 1 + pd_tol)) {
      stop(label, ": AR coefficients are not stationary; roots must lie outside the unit circle.",
           call. = FALSE)
    }
    invisible(TRUE)
  }

  as_autocor_spec <- function(x) {
    if (is.null(x)) return(list(structure = "NONE"))
    if (is.character(x)) return(list(structure = normalize_autocor(x)))
    if (!is.list(x)) {
      stop("Autocorrelation specification must be a character string or list.", call. = FALSE)
    }
    if (is.null(x$structure)) x$structure <- "NONE"
    x$structure <- normalize_autocor(x$structure)
    x
  }

  old_autocor_to_spec <- function(autocor, autocor_param) {
    spec <- list(structure = autocor)
    if (!is.null(autocor_param)) {
      if (is.list(autocor_param)) {
        spec <- c(spec, autocor_param)
      } else {
        spec$rho <- autocor_param
      }
    }
    as_autocor_spec(spec)
  }

  value_for_predictor <- function(x, j, nm, structure_j, param_name) {
    if (is.null(x)) return(NULL)
    if (is.matrix(x)) return(x)

    if (is.list(x)) {
      x_names <- names(x)
      if (!is.null(x_names) && nm %in% x_names) return(x[[nm]])
      if (length(x) == p && (is.null(x_names) || all(x_names == ""))) return(x[[j]])
      return(x)
    }

    if (!is.null(names(x)) && nm %in% names(x)) return(x[[nm]])

    if (param_name == "structure" && length(x) == p) return(x[[j]])

    scalar_structure <- structure_j %in% c("AR1", "CS")
    scalar_param <- param_name %in% c("rho", "phi")
    if (length(x) == p && scalar_structure && scalar_param) return(x[[j]])

    x
  }

  build_predictor_specs <- function(predictor_autocor) {
    if (is.null(predictor_autocor) && !is.null(phi)) {
      predictor_autocor <- list(structure = "AR1", rho = phi)
    }
    if (is.null(predictor_autocor)) {
      predictor_autocor <- list(structure = "NONE")
    }
    if (is.character(predictor_autocor)) {
      predictor_autocor <- list(structure = predictor_autocor)
    }
    if (!is.list(predictor_autocor)) {
      stop("predictor_autocor must be NULL, a character string, or a list.", call. = FALSE)
    }

    out <- vector("list", p)
    names(out) <- var_names

    if (!is.null(predictor_autocor$structure)) {
      for (j in seq_len(p)) {
        st_raw <- value_for_predictor(predictor_autocor$structure, j, var_names[j], "NONE", "structure")
        st <- normalize_autocor(st_raw)
        spec_j <- list(structure = st)
        for (par_name in setdiff(names(predictor_autocor), c("structure", "default", "by_predictor"))) {
          spec_j[[par_name]] <- value_for_predictor(
            predictor_autocor[[par_name]], j, var_names[j], st, par_name
          )
        }
        out[[j]] <- as_autocor_spec(spec_j)
      }
      return(out)
    }

    predictor_names <- names(predictor_autocor)
    if (is.null(predictor_names)) predictor_names <- rep("", length(predictor_autocor))
    unnamed_positions <- which(predictor_names == "" | is.na(predictor_names))
    default_spec <- predictor_autocor$default %||% list(structure = "NONE")

    for (j in seq_len(p)) {
      if (var_names[j] %in% predictor_names) {
        spec_j <- predictor_autocor[[var_names[j]]]
      } else if (as.character(j) %in% predictor_names) {
        spec_j <- predictor_autocor[[as.character(j)]]
      } else if (length(unnamed_positions) >= j) {
        spec_j <- predictor_autocor[[unnamed_positions[j]]]
      } else {
        spec_j <- default_spec
      }
      out[[j]] <- as_autocor_spec(spec_j)
    }

    out
  }

  autocor_matrix <- function(spec, label) {
    spec <- as_autocor_spec(spec)
    structure <- spec$structure

    if (n_time == 1L) {
      return(matrix(1, 1L, 1L))
    }

    lag_mat <- abs(outer(seq_len(n_time), seq_len(n_time), "-"))

    if (structure == "NONE") {
      R <- diag(n_time)
    } else if (structure == "AR1") {
      rho <- spec$rho %||% spec$phi
      if (is.null(rho) || length(rho) != 1L || !is.finite(rho)) {
        stop(label, ": AR1 requires one finite rho or phi.", call. = FALSE)
      }
      if (abs(rho) >= 1) {
        stop(label, ": AR1 requires abs(rho) < 1.", call. = FALSE)
      }
      R <- rho ^ lag_mat
    } else if (structure == "ARP") {
      ar <- spec$ar %||% spec$phi %||% spec$rho
      check_ar_stationarity(ar, label = label)
      acf_vec <- tryCatch(
        stats::ARMAacf(ar = as.numeric(ar), lag.max = n_time - 1L),
        error = function(e) {
          stop(label, ": invalid ARp coefficients: ", conditionMessage(e), call. = FALSE)
        }
      )
      R <- toeplitz(as.numeric(acf_vec))
    } else if (structure == "ARMA11") {
      ar <- spec$ar %||% spec$phi %||% spec$rho
      ma <- spec$ma %||% spec$theta
      if (is.null(ar) || length(ar) != 1L || !is.finite(ar) || abs(ar) >= 1) {
        stop(label, ": ARMA11 requires one finite AR coefficient with abs(ar) < 1.",
             call. = FALSE)
      }
      if (is.null(ma) || length(ma) != 1L || !is.finite(ma)) {
        stop(label, ": ARMA11 requires one finite MA coefficient named ma or theta.",
             call. = FALSE)
      }
      acf_vec <- tryCatch(
        stats::ARMAacf(ar = as.numeric(ar), ma = as.numeric(ma), lag.max = n_time - 1L),
        error = function(e) {
          stop(label, ": invalid ARMA11 coefficients: ", conditionMessage(e), call. = FALSE)
        }
      )
      R <- toeplitz(as.numeric(acf_vec))
    } else if (structure == "CS") {
      rho <- spec$rho %||% spec$phi
      if (is.null(rho) || length(rho) != 1L || !is.finite(rho)) {
        stop(label, ": CS requires one finite rho or phi.", call. = FALSE)
      }
      lower <- -1 / (n_time - 1L)
      if (rho <= lower || rho >= 1) {
        stop(label, ": CS requires -1/(n_time - 1) < rho < 1.", call. = FALSE)
      }
      R <- matrix(rho, n_time, n_time)
      diag(R) <- 1
    } else if (structure == "TOEP") {
      rho_vec <- spec$rho_vec %||% spec$rho %||% spec$phi
      if (is.null(rho_vec) || length(rho_vec) != n_time - 1L || any(!is.finite(rho_vec))) {
        stop(label, ": TOEP requires finite rho_vec, rho, or phi of length n_time - 1.",
             call. = FALSE)
      }
      if (any(abs(rho_vec) > 1)) {
        stop(label, ": TOEP correlations must lie in [-1, 1].", call. = FALSE)
      }
      R <- toeplitz(c(1, as.numeric(rho_vec)))
    } else if (structure %in% c("UN", "USER")) {
      R <- spec$R %||% spec$Sigma %||% spec$matrix %||% spec$corr
      if (is.null(R)) {
        stop(label, ": UN/USER requires a matrix named R, Sigma, matrix, or corr.",
             call. = FALSE)
      }
      R <- as.matrix(R)
    } else {
      stop(label, ": unsupported autocorrelation structure.", call. = FALSE)
    }

    as_positive_definite_corr(R, name = label)
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 3. Set simulator -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  predictor_specs <- build_predictor_specs(predictor_autocor)
  predictor_R <- vector("list", p)
  names(predictor_R) <- var_names
  for (j in seq_len(p)) {
    predictor_R[[j]] <- autocor_matrix(predictor_specs[[j]], label = paste0("predictor_autocor for ", var_names[j]))
  }
  has_predictor_autocor <- any(vapply(predictor_specs, function(x) x$structure != "NONE", logical(1L)))

  residual_spec <- if (is.null(residual_autocor)) {
    old_autocor_to_spec(autocor, autocor_param)
  } else {
    as_autocor_spec(residual_autocor)
  }
  R_eps <- autocor_matrix(residual_spec, label = "residual_autocor")
  Sigma_eps <- sigma_eps^2 * R_eps


  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 4. Generate predictors -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  joint_predictor_covariance <- NULL

  if (!has_predictor_autocor) {
    if (exact_predictor_cor) {
      if (n_obs <= p) {
        stop("exact_predictor_cor = TRUE requires n_id * n_time > length(beta).",
             call. = FALSE)
      }
      Z <- matrix(stats::rnorm(n_obs * p), nrow = n_obs, ncol = p)
      Z <- scale(Z, center = TRUE, scale = FALSE)
      S_Z <- crossprod(Z) / (n_obs - 1L)
      U_Z <- chol(S_Z)
      Z_white <- Z %*% solve(U_Z)
      X_std <- Z_white %*% chol(cor_matrix)
    } else {
      Z <- matrix(stats::rnorm(n_obs * p), nrow = n_obs, ncol = p)
      X_std <- Z %*% chol(cor_matrix)
    }

    X_mat <- sweep(X_std, 2L, predictor_sds, "*")
    X_mat <- sweep(X_mat, 2L, predictor_means, "+")
  } else {
    q <- p * n_time
    idx <- function(j) ((j - 1L) * n_time + 1L):(j * n_time)

    joint_predictor_covariance <- matrix(0, q, q)
    for (j in seq_len(p)) {
      for (k in seq_len(p)) {
        if (j == k) {
          block <- predictor_R[[j]]
        } else if (cross_lag_rule == "same_time_only") {
          block <- diag(n_time)
        } else if (cross_lag_rule == "zero") {
          block <- matrix(0, n_time, n_time)
        } else if (cross_lag_rule == "geometric") {
          block <- sign(predictor_R[[j]] * predictor_R[[k]]) *
            sqrt(abs(predictor_R[[j]] * predictor_R[[k]]))
        } else if (cross_lag_rule == "average") {
          block <- (predictor_R[[j]] + predictor_R[[k]]) / 2
        }

        joint_predictor_covariance[idx(j), idx(k)] <-
          cor_matrix[j, k] * predictor_sds[j] * predictor_sds[k] * block
      }
    }
    joint_predictor_covariance <- as_positive_definite_cov(
      joint_predictor_covariance,
      name = "joint predictor covariance matrix"
    )

    mean_x <- rep(predictor_means, each = n_time)

    if (exact_predictor_cor && n_id > q) {
      Z <- matrix(stats::rnorm(n_id * q), nrow = n_id, ncol = q)
      Z <- scale(Z, center = TRUE, scale = FALSE)
      S_Z <- crossprod(Z) / (n_id - 1L)
      U_Z <- chol(S_Z)
      Z_white <- Z %*% solve(U_Z)
      X_wide <- Z_white %*% chol(joint_predictor_covariance)
    } else {
      if (exact_predictor_cor && n_id <= q) {
        warning(
          "exact_predictor_cor = TRUE with predictor_autocor requires n_id > length(beta) * n_time. ",
          "Using population-target predictor covariance instead.",
          call. = FALSE
        )
      }
      X_wide <- matrix(stats::rnorm(n_id * q), nrow = n_id, ncol = q) %*%
        chol(joint_predictor_covariance)
    }
    X_wide <- sweep(X_wide, 2L, mean_x, "+")

    X_array <- array(NA_real_, dim = c(n_id, n_time, p))
    for (j in seq_len(p)) {
      X_array[, , j] <- X_wide[, idx(j), drop = FALSE]
    }

    X_mat <- matrix(NA_real_, nrow = n_obs, ncol = p)
    for (j in seq_len(p)) {
      X_mat[, j] <- as.vector(t(X_array[, , j]))
    }
  }

  colnames(X_mat) <- var_names


  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 5. Construct linear and nonlinear fixed-effect components -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  eta <- rep(intercept, n_obs) + as.numeric(X_mat %*% beta)

  if (!is.null(nonlinear_fun)) {
    if (!is.function(nonlinear_fun)) {
      stop("nonlinear_fun must be a function or NULL.", call. = FALSE)
    }
    nl <- nonlinear_fun(X_mat)
    if (!is.numeric(nl) || length(nl) != n_obs || any(!is.finite(nl))) {
      stop("nonlinear_fun must return a finite numeric vector of length n_id * n_time.",
           call. = FALSE)
    }
    eta <- eta + as.numeric(nl)
  } else if (nonlinear_type != "none") {
    if (nonlinear_type == "quadratic") {
      if (!is.numeric(beta_nl) || length(beta_nl) != p || any(!is.finite(beta_nl))) {
        stop("For nonlinear_type = 'quadratic', beta_nl must be a finite numeric vector of length length(beta).",
             call. = FALSE)
      }
      eta <- eta + as.numeric((X_mat^2) %*% beta_nl)
    }

    if (nonlinear_type == "interaction") {
      if (p < 2L) {
        stop("nonlinear_type = 'interaction' requires at least two predictors.", call. = FALSE)
      }
      pair_index <- utils::combn(seq_len(p), 2L)
      X_int <- X_mat[, pair_index[1L, ], drop = FALSE] *
        X_mat[, pair_index[2L, ], drop = FALSE]
      if (!is.numeric(beta_nl) || length(beta_nl) != ncol(X_int) || any(!is.finite(beta_nl))) {
        stop("For nonlinear_type = 'interaction', beta_nl must be a finite numeric vector of length choose(length(beta), 2).",
             call. = FALSE)
      }
      eta <- eta + as.numeric(X_int %*% beta_nl)
    }
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 6. Add subject-specific random intercepts -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  b_i <- if (sigma_b == 0) {
    rep(0, n_id)
  } else {
    stats::rnorm(n_id, mean = 0, sd = sigma_b)
  }
  eta <- eta + rep(b_i, each = n_time)

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 7. Generate latent residuals with selected residual autocorrelation -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  if (sigma_eps == 0) {
    latent_eps <- rep(0, n_obs)
  } else {
    U_eps <- chol(Sigma_eps)
    E <- matrix(stats::rnorm(n_id * n_time), nrow = n_id, ncol = n_time)
    E <- E %*% U_eps
    latent_eps <- as.vector(t(E))
  }
  eta <- eta + latent_eps

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 8. Generate observed outcome -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  if (family == "gaussian") {
    y <- eta
  } else if (family == "binomial") {
    prob <- stats::plogis(eta)
    y <- stats::rbinom(n = n_obs, size = binomial_size, prob = prob)
  } else if (family == "poisson") {
    lambda <- exp(eta)
    if (any(!is.finite(lambda))) {
      stop("Poisson rate overflow occurred. Reduce intercept, beta, sigma_b, or sigma_eps.",
           call. = FALSE)
    }
    y <- stats::rpois(n = n_obs, lambda = lambda)
  }

  dat <- data.frame(
    id = rep(seq_len(n_id), each = n_time),
    time = rep(seq_len(n_time), times = n_id),
    y = y,
    X_mat,
    check.names = FALSE
  )

  if (return_latent) {
    dat$eta <- eta
    dat$latent_eps <- latent_eps
  }

  rownames(dat) <- NULL

  if (return_components) {
    return(list(
      data = dat,
      predictor_specs = predictor_specs,
      predictor_autocor_matrices = predictor_R,
      residual_spec = residual_spec,
      residual_autocor_matrix = R_eps,
      joint_predictor_covariance = joint_predictor_covariance,
      eta = eta,
      latent_eps = latent_eps,
      random_intercepts = b_i,
      beta = beta,
      predictor_means = predictor_means,
      predictor_sds = predictor_sds,
      cross_lag_rule = cross_lag_rule
    ))
  }

  dat
}
