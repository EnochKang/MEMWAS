# Internal MEMWAS helpers. These functions are intentionally not exported.
# Internal native-C++ backend helpers for MEMWAS. These functions intentionally
# use base R's .Call interface only; no Rcpp, RcppArmadillo, or external package
# dependency is required.

.memwas_normalize_engine <- function(engine = c("R", "cpp")) {
  if (is.null(engine)) engine <- "R"
  engine <- tolower(as.character(engine)[1L])
  aliases <- c(r = "R", base = "R", base_r = "R", rbase = "R",
               cpp = "cpp", `c++` = "cpp", cplusplus = "cpp", cxx = "cpp", native = "cpp")
  if (!engine %in% names(aliases)) {
    stop("`engine` must be either 'R' or 'cpp'.", call. = FALSE)
  }
  unname(aliases[[engine]])
}

.memwas_is_cpp <- function(engine) identical(.memwas_normalize_engine(engine), "cpp")

.memwas_safe_exp <- function(x, lo = -30, hi = 30, engine = "R") {
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_safe_exp", as.numeric(x), as.numeric(lo), as.numeric(hi), PACKAGE = "MEMWAS"))
  }
  exp(pmin(pmax(x, lo), hi))
}

.memwas_diag_vec <- function(v, engine = "R") {
  v <- as.numeric(v)
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_diag_vec", v, PACKAGE = "MEMWAS"))
  }
  if (length(v) == 0L) return(matrix(0, 0L, 0L))
  M <- matrix(0, length(v), length(v))
  diag(M) <- v
  M
}

.memwas_bounded_logistic <- function(x, low, high, engine = "R") {
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_bounded_logistic", as.numeric(x), as.numeric(low), as.numeric(high), PACKAGE = "MEMWAS"))
  }
  z <- 1 / (1 + exp(-pmin(pmax(x, -30), 30)))
  low + (high - low) * z
}

.memwas_soft_threshold <- function(z, gamma, engine = "R") {
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_soft_threshold", as.numeric(z), as.numeric(gamma), PACKAGE = "MEMWAS"))
  }
  sign(z) * pmax(abs(z) - gamma, 0)
}

.memwas_penalty_factor <- function(xnames, engine = "R") {
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_penalty_factor", as.character(xnames), PACKAGE = "MEMWAS"))
  }
  pf <- rep(1, length(xnames))
  pf[xnames == "(Intercept)"] <- 0
  pf
}

.memwas_fixed_effect_penalty <- function(beta, pf, lambda1, lambda2, engine = "R") {
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_fixed_effect_penalty", as.numeric(beta), as.numeric(pf),
                 as.numeric(lambda1), as.numeric(lambda2), PACKAGE = "MEMWAS"))
  }
  lambda1 * sum(abs(beta) * pf) + 0.5 * lambda2 * sum((beta^2) * pf)
}

.memwas_safe_solve <- function(A, b = NULL, eps = 1e-8, engine = "R") {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("Internal error: solve target is not square.", call. = FALSE)
  if (is.null(b)) b <- diag(nrow(A))
  if (.memwas_is_cpp(engine)) {
    B <- as.matrix(b)
    storage.mode(A) <- "double"
    storage.mode(B) <- "double"
    return(.Call("memwas_safe_solve", A, B, as.numeric(eps), PACKAGE = "MEMWAS"))
  }
  A <- (A + t(A)) / 2
  jitters <- c(0, eps, eps * 10, eps * 1000, eps * 1e5)
  for (j in jitters) {
    out <- try(solve(A + diag(j, nrow(A)), b), silent = TRUE)
    if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  }
  out <- try(qr.solve(A + diag(eps * 1e5, nrow(A)), b), silent = TRUE)
  if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  stop("Matrix solve failed; design or covariance matrix is numerically singular.", call. = FALSE)
}

.memwas_safe_chol <- function(A, eps = 1e-8, engine = "R") {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("Internal error: Cholesky target is not square.", call. = FALSE)
  if (.memwas_is_cpp(engine)) {
    storage.mode(A) <- "double"
    return(.Call("memwas_safe_chol", A, as.numeric(eps), PACKAGE = "MEMWAS"))
  }
  A <- (A + t(A)) / 2
  if (nrow(A) == 1L) {
    return(matrix(sqrt(max(A[1L, 1L], eps)), 1L, 1L))
  }
  jitters <- c(0, eps, eps * 10, eps * 1000, eps * 1e5)
  for (j in jitters) {
    out <- try(chol(A + diag(j, nrow(A))), silent = TRUE)
    if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  }
  eg <- eigen(A, symmetric = TRUE)
  vals <- pmax(eg$values, eps)
  A2 <- eg$vectors %*% diag(vals, nrow = length(vals)) %*% t(eg$vectors)
  A2 <- (A2 + t(A2)) / 2
  chol(A2 + diag(eps, nrow(A2)))
}

.memwas_coordinate_descent_enet <- function(X, y, lambda1, lambda2, pf,
                                            maxit = 1000L, tol = 1e-7,
                                            engine = "R") {
  X <- as.matrix(X)
  y <- as.numeric(y)
  pf <- as.numeric(pf)
  p <- ncol(X)
  if (p == 0L) return(numeric(0L))
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_coordinate_descent_enet", X, y, as.numeric(lambda1),
                 as.numeric(lambda2), pf, as.integer(maxit), as.numeric(tol),
                 PACKAGE = "MEMWAS"))
  }
  H <- crossprod(X) + .memwas_diag_vec(lambda2 * pf + 1e-8, engine = "R")
  beta <- as.vector(.memwas_safe_solve(H, crossprod(X, y), engine = "R"))
  x2 <- colSums(X * X)
  r <- as.vector(y - X %*% beta)
  for (iter in seq_len(maxit)) {
    beta_old <- beta
    for (j in seq_len(p)) {
      r <- r + X[, j] * beta[j]
      zj <- sum(X[, j] * r)
      denom <- x2[j] + lambda2 * pf[j] + 1e-12
      if (pf[j] == 0) {
        beta[j] <- zj / denom
      } else {
        beta[j] <- .memwas_soft_threshold(zj, lambda1 * pf[j], engine = "R") / denom
      }
      r <- r - X[, j] * beta[j]
    }
    if (max(abs(beta - beta_old)) < tol) break
  }
  beta
}

.memwas_rcs_basis <- function(x, knots, engine = "R") {
  x <- as.numeric(x)
  knots <- sort(unique(as.numeric(knots)))
  if (.memwas_is_cpp(engine)) {
    B <- .Call("memwas_rcs_basis", x, knots, PACKAGE = "MEMWAS")
    colnames(B) <- paste0("spline", seq_len(ncol(B)))
    return(B)
  }
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

.memwas_lag_matrix <- function(ti, continuous = FALSE, engine = "R") {
  n <- length(ti)
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_lag_matrix", suppressWarnings(as.numeric(ti)), as.logical(continuous), PACKAGE = "MEMWAS"))
  }
  if (n <= 1L) return(matrix(0, n, n))
  if (continuous) {
    x <- if (inherits(ti, "POSIXt") || inherits(ti, "Date")) as.numeric(ti) else suppressWarnings(as.numeric(ti))
    if (all(is.finite(x))) {
      ux <- sort(unique(x))
      dif <- diff(ux)
      scale <- if (length(dif) && any(dif > 0)) min(dif[dif > 0]) else 1
      return(abs(outer(x, x, "-")) / scale)
    }
  }
  abs(outer(seq_len(n), seq_len(n), "-"))
}

.memwas_pacf_to_ar <- function(pacf, engine = "R") {
  pacf <- as.numeric(pacf)
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_pacf_to_ar", pacf, PACKAGE = "MEMWAS"))
  }
  p <- length(pacf)
  if (p == 0L) return(numeric(0L))
  phi <- matrix(0, p, p)
  for (k in seq_len(p)) {
    phi[k, k] <- pacf[k]
    if (k > 1L) {
      for (j in seq_len(k - 1L)) {
        phi[k, j] <- phi[k - 1L, j] - pacf[k] * phi[k - 1L, k - j]
      }
    }
  }
  as.numeric(phi[p, seq_len(p)])
}

.memwas_ar_acf <- function(phi, max_lag, engine = "R") {
  phi <- as.numeric(phi)
  max_lag <- as.integer(max_lag)
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_ar_acf", phi, max_lag, PACKAGE = "MEMWAS"))
  }
  p <- length(phi)
  if (max_lag == 0L) return(1)
  A <- diag(1, p)
  cvec <- numeric(p)
  for (k in seq_len(p)) {
    for (j in seq_len(p)) {
      lag <- abs(k - j)
      if (lag == 0L) cvec[k] <- cvec[k] + phi[j] else A[k, lag] <- A[k, lag] - phi[j]
    }
  }
  rho <- numeric(max(max_lag, p) + 1L)
  rho[1L] <- 1
  rho[2L:(p + 1L)] <- as.numeric(.memwas_safe_solve(A, cvec, engine = "R"))
  if (max_lag > p) {
    for (k in (p + 1L):max_lag) rho[k + 1L] <- sum(phi * rho[k - seq_len(p) + 1L])
  }
  pmax(pmin(rho[seq_len(max_lag + 1L)], 0.999), -0.999)
}

.memwas_kernel_matrix <- function(A, B, ell, engine = "R") {
  A <- as.matrix(A)
  B <- as.matrix(B)
  ell <- as.numeric(ell)
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_kernel_matrix", A, B, ell, PACKAGE = "MEMWAS"))
  }
  out <- matrix(0, nrow(A), nrow(B))
  for (j in seq_len(ncol(A))) {
    out <- out + (outer(A[, j], B[, j], "-") / ell[j])^2
  }
  exp(-0.5 * out)
}

.memwas_metric_value <- function(y, pred, metric, epsilon = .Machine$double.eps,
                                 fail_value = Inf, engine = "R") {
  metric <- toupper(as.character(metric)[1L])
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_metric_value", as.numeric(y), as.numeric(pred), metric,
                 as.numeric(epsilon), as.numeric(fail_value), PACKAGE = "MEMWAS"))
  }
  ok <- is.finite(y) & is.finite(pred)
  y <- y[ok]
  pred <- pred[ok]
  if (length(y) == 0L) return(fail_value)
  e <- y - pred
  if (metric == "MAE") return(mean(abs(e)))
  if (metric == "MSE") return(mean(e^2))
  if (metric == "RMSE") return(sqrt(mean(e^2)))
  if (metric == "MAPE") return(100 * mean(abs(e) / pmax(abs(y), epsilon)))
  if (metric == "SMAPE") return(100 * mean(2 * abs(e) / pmax(abs(y) + abs(pred), epsilon)))
  stop("Internal error: unsupported metric.", call. = FALSE)
}

.memwas_stability_value <- function(x, metric, engine = "R") {
  metric <- toupper(as.character(metric)[1L])
  if (metric == "VAR") metric <- "VARIANCE"
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_stability_value", as.numeric(x), metric, PACKAGE = "MEMWAS"))
  }
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(0)
  if (metric == "IQR") return(as.numeric(stats::IQR(x, na.rm = TRUE)))
  if (metric == "SD") return(as.numeric(stats::sd(x, na.rm = TRUE)))
  if (metric == "VARIANCE") return(as.numeric(stats::var(x, na.rm = TRUE)))
  stop("Internal error: unsupported stability metric.", call. = FALSE)
}

.memwas_make_fold_assignment <- function(id, K, engine = "R") {
  K <- as.integer(K)
  if (.memwas_is_cpp(engine)) {
    return(.Call("memwas_make_fold_assignment", as.character(id), K, PACKAGE = "MEMWAS"))
  }
  counts <- table(as.character(id))
  n_id <- length(counts)
  if (K > n_id) {
    stop("`K` cannot exceed the number of unique non-missing subject IDs.", call. = FALSE)
  }
  ids_random <- sample(names(counts), n_id, replace = FALSE)
  fold_load <- rep(0, K)
  assigned <- integer(n_id)
  names(assigned) <- ids_random
  for (sid in ids_random) {
    fold <- which.min(fold_load)
    assigned[sid] <- fold
    fold_load[fold] <- fold_load[fold] + as.integer(counts[[sid]])
  }
  as.integer(assigned[as.character(id)])
}
