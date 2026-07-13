
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <climits>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef R_NO_REMAP
# define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

static double as_scalar_real(SEXP x, double default_value = 0.0) {
  if (Rf_length(x) < 1) return default_value;
  switch (TYPEOF(x)) {
  case REALSXP: return REAL(x)[0];
  case INTSXP:
  case LGLSXP: return INTEGER(x)[0] == NA_INTEGER ? NA_REAL : static_cast<double>(INTEGER(x)[0]);
  default: return Rf_asReal(x);
  }
}

static int as_scalar_int(SEXP x, int default_value = 0) {
  if (Rf_length(x) < 1) return default_value;
  switch (TYPEOF(x)) {
  case INTSXP:
  case LGLSXP: return INTEGER(x)[0];
  case REALSXP: return static_cast<int>(REAL(x)[0]);
  default: return Rf_asInteger(x);
  }
}

static inline R_xlen_t idx(int i, int j, int nrow) {
  return static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(j) * nrow;
}

static double clamp_double(double x, double lo, double hi) {
  if (!R_FINITE(x)) return x;
  return std::min(std::max(x, lo), hi);
}

static double soft_threshold_scalar(double z, double gamma) {
  if (z > gamma) return z - gamma;
  if (z < -gamma) return z + gamma;
  return 0.0;
}

static double max_abs_value(const std::vector<double>& x) {
  double out = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (!R_FINITE(x[i])) return std::numeric_limits<double>::infinity();
    out = std::max(out, std::fabs(x[i]));
  }
  return out;
}

static bool all_finite(const std::vector<double>& x) {
  for (size_t i = 0; i < x.size(); ++i) {
    if (!R_FINITE(x[i])) return false;
  }
  return true;
}

static void symmetrize_square(std::vector<double>& A, int n) {
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      double aij = A[idx(i, j, n)];
      double aji = A[idx(j, i, n)];
      double s = 0.5 * (aij + aji);
      A[idx(i, j, n)] = s;
      A[idx(j, i, n)] = s;
    }
  }
}

static void swap_rows(std::vector<double>& M, int nrow, int ncol, int r1, int r2) {
  if (r1 == r2) return;
  for (int j = 0; j < ncol; ++j) {
    std::swap(M[idx(r1, j, nrow)], M[idx(r2, j, nrow)]);
  }
}

static bool gaussian_elimination_solve(std::vector<double> A,
                                       std::vector<double>& B,
                                       int n,
                                       int nrhs,
                                       double pivot_tol) {
  if (n == 0) return true;
  if (static_cast<int>(B.size()) != n * nrhs) return false;

  for (int k = 0; k < n; ++k) {
    int pivot = k;
    double pivot_abs = std::fabs(A[idx(k, k, n)]);
    for (int i = k + 1; i < n; ++i) {
      double cand = std::fabs(A[idx(i, k, n)]);
      if (cand > pivot_abs) {
        pivot = i;
        pivot_abs = cand;
      }
    }

    if (!R_FINITE(pivot_abs) || pivot_abs <= pivot_tol) return false;

    if (pivot != k) {
      swap_rows(A, n, n, pivot, k);
      swap_rows(B, n, nrhs, pivot, k);
    }

    double akk = A[idx(k, k, n)];
    if (!R_FINITE(akk) || std::fabs(akk) <= pivot_tol) return false;

    for (int i = k + 1; i < n; ++i) {
      double factor = A[idx(i, k, n)] / akk;
      if (!R_FINITE(factor)) return false;
      A[idx(i, k, n)] = 0.0;
      if (factor == 0.0) continue;
      for (int j = k + 1; j < n; ++j) {
        A[idx(i, j, n)] -= factor * A[idx(k, j, n)];
      }
      for (int r = 0; r < nrhs; ++r) {
        B[idx(i, r, n)] -= factor * B[idx(k, r, n)];
      }
    }
  }

  for (int r = 0; r < nrhs; ++r) {
    for (int i = n - 1; i >= 0; --i) {
      double sum = B[idx(i, r, n)];
      for (int j = i + 1; j < n; ++j) {
        sum -= A[idx(i, j, n)] * B[idx(j, r, n)];
      }
      double aii = A[idx(i, i, n)];
      if (!R_FINITE(aii) || std::fabs(aii) <= pivot_tol || !R_FINITE(sum)) return false;
      B[idx(i, r, n)] = sum / aii;
    }
  }

  return all_finite(B);
}

static bool cholesky_upper(const std::vector<double>& A,
                           std::vector<double>& U,
                           int n,
                           double diag_tol) {
  U.assign(static_cast<size_t>(n) * n, 0.0);
  if (n == 0) return true;

  for (int j = 0; j < n; ++j) {
    double d = A[idx(j, j, n)];
    for (int k = 0; k < j; ++k) {
      double ukj = U[idx(k, j, n)];
      d -= ukj * ukj;
    }
    if (!R_FINITE(d) || d <= diag_tol) return false;
    U[idx(j, j, n)] = std::sqrt(d);

    for (int col = j + 1; col < n; ++col) {
      double s = A[idx(j, col, n)];
      for (int k = 0; k < j; ++k) {
        s -= U[idx(k, j, n)] * U[idx(k, col, n)];
      }
      if (!R_FINITE(s)) return false;
      U[idx(j, col, n)] = s / U[idx(j, j, n)];
    }
  }

  return all_finite(U);
}

static bool solve_with_jitter(const std::vector<double>& A0,
                              const std::vector<double>& B0,
                              int n,
                              int nrhs,
                              double eps,
                              bool symmetrize,
                              std::vector<double>& solution) {
  solution = B0;
  if (n == 0) return true;

  std::vector<double> baseA = A0;
  if (symmetrize) symmetrize_square(baseA, n);

  double scale = max_abs_value(baseA);
  if (!R_FINITE(scale) || scale < 1.0) scale = 1.0;
  double base_jitter = std::max(eps, 1e-12) * scale;
  double pivot_tol = std::max(1e-14, 100.0 * std::numeric_limits<double>::epsilon() * scale);
  const double multipliers[] = {0.0, 1.0, 10.0, 100.0, 1000.0, 1e5, 1e7, 1e9, 1e11};

  for (size_t trial = 0; trial < sizeof(multipliers) / sizeof(multipliers[0]); ++trial) {
    std::vector<double> A = baseA;
    std::vector<double> B = B0;
    double jitter = base_jitter * multipliers[trial];
    if (jitter > 0.0) {
      for (int i = 0; i < n; ++i) A[idx(i, i, n)] += jitter;
    }
    if (gaussian_elimination_solve(A, B, n, nrhs, pivot_tol)) {
      solution.swap(B);
      return true;
    }
  }

  return false;
}

static bool chol_with_jitter(const std::vector<double>& A0,
                             int n,
                             double eps,
                             std::vector<double>& U) {
  U.clear();
  if (n == 0) return true;

  std::vector<double> baseA = A0;
  symmetrize_square(baseA, n);

  double scale = max_abs_value(baseA);
  if (!R_FINITE(scale) || scale < 1.0) scale = 1.0;
  double base_jitter = std::max(eps, 1e-12) * scale;
  double diag_tol = std::max(1e-14, 100.0 * std::numeric_limits<double>::epsilon() * scale);
  const double multipliers[] = {0.0, 1.0, 10.0, 100.0, 1000.0, 1e5, 1e7, 1e9, 1e11};

  for (size_t trial = 0; trial < sizeof(multipliers) / sizeof(multipliers[0]); ++trial) {
    std::vector<double> A = baseA;
    double jitter = base_jitter * multipliers[trial];
    if (jitter > 0.0) {
      for (int i = 0; i < n; ++i) A[idx(i, i, n)] += jitter;
    }
    if (cholesky_upper(A, U, n, diag_tol)) {
      for (int j = 0; j < n; ++j) {
        for (int i = j + 1; i < n; ++i) U[idx(i, j, n)] = 0.0;
      }
      return true;
    }
  }

  return false;
}

extern "C" {
  SEXP memwas_safe_exp(SEXP xSEXP, SEXP loSEXP, SEXP hiSEXP) {
    SEXP x = PROTECT(Rf_coerceVector(xSEXP, REALSXP));
    R_xlen_t n = Rf_xlength(x);
    double lo = as_scalar_real(loSEXP, -30.0);
    double hi = as_scalar_real(hiSEXP, 30.0);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    const double* px = REAL(x);
    double* po = REAL(out);
    for (R_xlen_t i = 0; i < n; ++i) po[i] = std::exp(clamp_double(px[i], lo, hi));
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_diag_vec(SEXP vSEXP) {
    SEXP v = PROTECT(Rf_coerceVector(vSEXP, REALSXP));
    R_xlen_t n0 = Rf_xlength(v);
    if (n0 > INT_MAX) Rf_error("Vector is too long to form a diagonal matrix.");
    int n = static_cast<int>(n0);
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    std::fill(REAL(out), REAL(out) + static_cast<R_xlen_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) REAL(out)[idx(i, i, n)] = REAL(v)[i];
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_bounded_logistic(SEXP xSEXP, SEXP lowSEXP, SEXP highSEXP) {
    SEXP x = PROTECT(Rf_coerceVector(xSEXP, REALSXP));
    double low = as_scalar_real(lowSEXP, 0.0);
    double high = as_scalar_real(highSEXP, 1.0);
    R_xlen_t n = Rf_xlength(x);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
      double z = 1.0 / (1.0 + std::exp(-clamp_double(REAL(x)[i], -30.0, 30.0)));
      REAL(out)[i] = low + (high - low) * z;
    }
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_soft_threshold(SEXP zSEXP, SEXP gammaSEXP) {
    SEXP z = PROTECT(Rf_coerceVector(zSEXP, REALSXP));
    double gamma = as_scalar_real(gammaSEXP, 0.0);
    R_xlen_t n = Rf_xlength(z);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) REAL(out)[i] = soft_threshold_scalar(REAL(z)[i], gamma);
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_penalty_factor(SEXP namesSEXP) {
    SEXP nm = PROTECT(Rf_coerceVector(namesSEXP, STRSXP));
    R_xlen_t n = Rf_xlength(nm);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
      const char* s = CHAR(STRING_ELT(nm, i));
      REAL(out)[i] = (std::string(s) == "(Intercept)") ? 0.0 : 1.0;
    }
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_fixed_effect_penalty(SEXP betaSEXP, SEXP pfSEXP, SEXP l1SEXP, SEXP l2SEXP) {
    SEXP beta = PROTECT(Rf_coerceVector(betaSEXP, REALSXP));
    SEXP pf = PROTECT(Rf_coerceVector(pfSEXP, REALSXP));
    R_xlen_t n = Rf_xlength(beta);
    if (Rf_xlength(pf) != n) Rf_error("Penalty factor length must equal coefficient length.");
    double l1 = as_scalar_real(l1SEXP, 0.0);
    double l2 = as_scalar_real(l2SEXP, 0.0);
    double val = 0.0;
    for (R_xlen_t i = 0; i < n; ++i) {
      double b = REAL(beta)[i];
      double f = REAL(pf)[i];
      val += l1 * std::fabs(b) * f + 0.5 * l2 * b * b * f;
    }
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(out)[0] = val;
    UNPROTECT(3);
    return out;
  }
}

extern "C" {
  SEXP memwas_safe_solve(SEXP ASEXP, SEXP BSEXP, SEXP epsSEXP) {
    if (!Rf_isMatrix(ASEXP)) Rf_error("A must be a matrix.");

    SEXP dimA = Rf_getAttrib(ASEXP, R_DimSymbol);
    int n = INTEGER(dimA)[0];
    int ncolA = INTEGER(dimA)[1];
    if (n != ncolA) Rf_error("A must be square.");

    SEXP Ain = PROTECT(Rf_coerceVector(ASEXP, REALSXP));

    int nb = 0;
    int nrhs = 0;
    SEXP Bin = R_NilValue;
    bool b_is_matrix = Rf_isMatrix(BSEXP);
    if (b_is_matrix) {
      SEXP dimB = Rf_getAttrib(BSEXP, R_DimSymbol);
      nb = INTEGER(dimB)[0];
      nrhs = INTEGER(dimB)[1];
      if (nb != n) Rf_error("Non-conformable matrix right-hand side.");
      Bin = PROTECT(Rf_coerceVector(BSEXP, REALSXP));
    } else {
      Bin = PROTECT(Rf_coerceVector(BSEXP, REALSXP));
      if (Rf_xlength(Bin) != n) Rf_error("Non-conformable right-hand side.");
      nb = n;
      nrhs = 1;
    }

    std::vector<double> A(static_cast<size_t>(n) * n);
    std::copy(REAL(Ain), REAL(Ain) + static_cast<R_xlen_t>(n) * n, A.begin());

    std::vector<double> B(static_cast<size_t>(nb) * nrhs);
    std::copy(REAL(Bin), REAL(Bin) + static_cast<R_xlen_t>(nb) * nrhs, B.begin());

    double eps = as_scalar_real(epsSEXP, 1e-8);
    std::vector<double> solution;
    if (!solve_with_jitter(A, B, n, nrhs, eps, true, solution)) {
      UNPROTECT(2);
      Rf_error("C++ safe solve failed; matrix is numerically singular.");
    }

    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, nrhs));
    std::copy(solution.begin(), solution.end(), REAL(out));
    UNPROTECT(3);
    return out;
  }
}

extern "C" {
  SEXP memwas_safe_chol(SEXP ASEXP, SEXP epsSEXP) {
    if (!Rf_isMatrix(ASEXP)) Rf_error("A must be a matrix.");
    SEXP dimA = Rf_getAttrib(ASEXP, R_DimSymbol);
    int n = INTEGER(dimA)[0];
    int ncolA = INTEGER(dimA)[1];
    if (n != ncolA) Rf_error("A must be square.");

    SEXP Ain = PROTECT(Rf_coerceVector(ASEXP, REALSXP));
    double eps = as_scalar_real(epsSEXP, 1e-8);

    if (n == 1) {
      SEXP out1 = PROTECT(Rf_allocMatrix(REALSXP, 1, 1));
      double a11 = REAL(Ain)[0];
      REAL(out1)[0] = std::sqrt(std::max(a11, eps));
      UNPROTECT(2);
      return out1;
    }

    std::vector<double> A(static_cast<size_t>(n) * n);
    std::copy(REAL(Ain), REAL(Ain) + static_cast<R_xlen_t>(n) * n, A.begin());
    std::vector<double> U;
    if (!chol_with_jitter(A, n, eps, U)) {
      UNPROTECT(1);
      Rf_error("C++ safe Cholesky failed; matrix is not positive definite after jitter.");
    }

    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    std::copy(U.begin(), U.end(), REAL(out));
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_coordinate_descent_enet(SEXP XSEXP, SEXP ySEXP, SEXP l1SEXP,
                                      SEXP l2SEXP, SEXP pfSEXP,
                                      SEXP maxitSEXP, SEXP tolSEXP) {
    if (!Rf_isMatrix(XSEXP)) Rf_error("X must be a matrix.");
    SEXP X = PROTECT(Rf_coerceVector(XSEXP, REALSXP));
    SEXP y = PROTECT(Rf_coerceVector(ySEXP, REALSXP));
    SEXP pf = PROTECT(Rf_coerceVector(pfSEXP, REALSXP));
    SEXP dimX = Rf_getAttrib(XSEXP, R_DimSymbol);
    int n = INTEGER(dimX)[0];
    int p = INTEGER(dimX)[1];
    if (Rf_xlength(y) != n || Rf_xlength(pf) != p) Rf_error("Non-conformable coordinate descent inputs.");
    double l1 = as_scalar_real(l1SEXP, 0.0);
    double l2 = as_scalar_real(l2SEXP, 0.0);
    int maxit = as_scalar_int(maxitSEXP, 1000);
    double tol = as_scalar_real(tolSEXP, 1e-7);
    SEXP beta_result = PROTECT(Rf_allocVector(REALSXP, p));
    std::vector<double> beta(p, 0.0), beta_old(p, 0.0), residual(n), x2(p, 0.0);
    for (int i = 0; i < n; ++i) residual[i] = REAL(y)[i];
    for (int j = 0; j < p; ++j) {
      double ss = 0.0;
      for (int i = 0; i < n; ++i) {
        double xij = REAL(X)[idx(i, j, n)];
        ss += xij * xij;
      }
      x2[j] = ss;
    }
    for (int iter = 0; iter < maxit; ++iter) {
      beta_old = beta;
      for (int j = 0; j < p; ++j) {
        for (int i = 0; i < n; ++i) residual[i] += REAL(X)[idx(i, j, n)] * beta[j];
        double z = 0.0;
        for (int i = 0; i < n; ++i) z += REAL(X)[idx(i, j, n)] * residual[i];
        double pfj = REAL(pf)[j];
        double denom = x2[j] + l2 * pfj + 1e-12;
        beta[j] = (pfj == 0.0) ? z / denom : soft_threshold_scalar(z, l1 * pfj) / denom;
        for (int i = 0; i < n; ++i) residual[i] -= REAL(X)[idx(i, j, n)] * beta[j];
      }
      double maxdiff = 0.0;
      for (int j = 0; j < p; ++j) maxdiff = std::max(maxdiff, std::fabs(beta[j] - beta_old[j]));
      if (maxdiff < tol) break;
    }
    for (int j = 0; j < p; ++j) REAL(beta_result)[j] = beta[j];
    UNPROTECT(4);
    return beta_result;
  }
}

extern "C" {
  SEXP memwas_rcs_basis(SEXP xSEXP, SEXP knotsSEXP) {
    SEXP x = PROTECT(Rf_coerceVector(xSEXP, REALSXP));
    SEXP knots0 = PROTECT(Rf_coerceVector(knotsSEXP, REALSXP));
    std::vector<double> knots(REAL(knots0), REAL(knots0) + Rf_xlength(knots0));
    std::sort(knots.begin(), knots.end());
    knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
    int K = static_cast<int>(knots.size());
    if (K < 3) Rf_error("At least three unique knots are required for restricted cubic splines.");
    R_xlen_t n = Rf_xlength(x);
    int nb = K - 2;
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n), nb));
    auto tp = [](double z, double k) { double v = z - k; return v > 0.0 ? v * v * v : 0.0; };
    double denom = std::max(knots[K - 1] - knots[K - 2], 1e-12);
    double scale = std::pow(std::max(knots[K - 1] - knots[0], 1e-12), 3.0);
    for (int j = 0; j < nb; ++j) {
      double kj = knots[j];
      for (R_xlen_t i = 0; i < n; ++i) {
        double xi = REAL(x)[i];
        double val = tp(xi, kj) - tp(xi, knots[K - 2]) * ((knots[K - 1] - kj) / denom) +
          tp(xi, knots[K - 1]) * ((knots[K - 2] - kj) / denom);
        REAL(out)[idx(static_cast<int>(i), j, static_cast<int>(n))] = val / scale;
      }
    }
    UNPROTECT(3);
    return out;
  }
}

extern "C" {
  SEXP memwas_lag_matrix(SEXP tSEXP, SEXP continuousSEXP) {
    SEXP t = PROTECT(Rf_coerceVector(tSEXP, REALSXP));
    int n = static_cast<int>(Rf_xlength(t));
    bool continuous = LOGICAL(continuousSEXP)[0] == TRUE;
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    bool all_finite = true;
    for (int i = 0; i < n; ++i) {
      if (!R_FINITE(REAL(t)[i])) {
        all_finite = false;
        break;
      }
    }
    double scale = 1.0;
    if (continuous && all_finite && n > 1) {
      std::vector<double> u(REAL(t), REAL(t) + n);
      std::sort(u.begin(), u.end());
      u.erase(std::unique(u.begin(), u.end()), u.end());
      double minpos = std::numeric_limits<double>::infinity();
      for (size_t i = 1; i < u.size(); ++i) {
        double d = u[i] - u[i - 1];
        if (d > 0.0 && d < minpos) minpos = d;
      }
      if (R_FINITE(minpos)) scale = minpos;
    }
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        double val = 0.0;
        if (continuous && all_finite) {
          val = std::fabs(REAL(t)[i] - REAL(t)[j]) / scale;
        } else {
          val = std::fabs(static_cast<double>((i + 1) - (j + 1)));
        }
        REAL(out)[idx(i, j, n)] = val;
      }
    }
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_pacf_to_ar(SEXP pacfSEXP) {
    SEXP pacf = PROTECT(Rf_coerceVector(pacfSEXP, REALSXP));
    int p = static_cast<int>(Rf_xlength(pacf));
    SEXP out = PROTECT(Rf_allocVector(REALSXP, p));
    if (p == 0) {
      UNPROTECT(2);
      return out;
    }
    std::vector<double> phi(static_cast<size_t>(p) * p, 0.0);
    for (int k = 0; k < p; ++k) {
      phi[idx(k, k, p)] = REAL(pacf)[k];
      if (k > 0) {
        for (int j = 0; j < k; ++j) {
          phi[idx(k, j, p)] = phi[idx(k - 1, j, p)] - REAL(pacf)[k] * phi[idx(k - 1, k - 1 - j, p)];
        }
      }
    }
    for (int j = 0; j < p; ++j) REAL(out)[j] = phi[idx(p - 1, j, p)];
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_ar_acf(SEXP phiSEXP, SEXP maxLagSEXP) {
    SEXP phiS = PROTECT(Rf_coerceVector(phiSEXP, REALSXP));
    int p = static_cast<int>(Rf_xlength(phiS));
    int max_lag = as_scalar_int(maxLagSEXP, 0);
    if (max_lag < 0) max_lag = 0;
    SEXP out = PROTECT(Rf_allocVector(REALSXP, max_lag + 1));

    if (max_lag == 0) {
      REAL(out)[0] = 1.0;
      UNPROTECT(2);
      return out;
    }
    if (p == 0) {
      REAL(out)[0] = 1.0;
      for (int i = 1; i <= max_lag; ++i) REAL(out)[i] = 0.0;
      UNPROTECT(2);
      return out;
    }

    std::vector<double> A(static_cast<size_t>(p) * p, 0.0), cvec(p, 0.0);
    for (int i = 0; i < p; ++i) A[idx(i, i, p)] = 1.0;
    for (int k = 0; k < p; ++k) {
      for (int j = 0; j < p; ++j) {
        int lag = std::abs((k + 1) - (j + 1));
        if (lag == 0) cvec[k] += REAL(phiS)[j];
        else A[idx(k, lag - 1, p)] -= REAL(phiS)[j];
      }
    }

    std::vector<double> rhs = cvec;
    bool ok = solve_with_jitter(A, rhs, p, 1, 1e-10, false, rhs);
    std::vector<double> rho(static_cast<size_t>(std::max(max_lag, p)) + 1, 0.0);
    rho[0] = 1.0;
    if (ok) {
      for (int i = 0; i < p; ++i) rho[i + 1] = rhs[i];
    }
    if (max_lag > p) {
      for (int k = p + 1; k <= max_lag; ++k) {
        double s = 0.0;
        for (int j = 1; j <= p; ++j) s += REAL(phiS)[j - 1] * rho[k - j];
        rho[k] = s;
      }
    }
    for (int i = 0; i <= max_lag; ++i) REAL(out)[i] = std::min(std::max(rho[i], -0.999), 0.999);
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_kernel_matrix(SEXP ASEXP, SEXP BSEXP, SEXP ellSEXP) {
    if (!Rf_isMatrix(ASEXP) || !Rf_isMatrix(BSEXP)) Rf_error("A and B must be matrices.");
    SEXP A = PROTECT(Rf_coerceVector(ASEXP, REALSXP));
    SEXP B = PROTECT(Rf_coerceVector(BSEXP, REALSXP));
    SEXP ell = PROTECT(Rf_coerceVector(ellSEXP, REALSXP));
    SEXP dimA = Rf_getAttrib(ASEXP, R_DimSymbol);
    SEXP dimB = Rf_getAttrib(BSEXP, R_DimSymbol);
    int nA = INTEGER(dimA)[0], dA = INTEGER(dimA)[1];
    int nB = INTEGER(dimB)[0], dB = INTEGER(dimB)[1];
    if (dA != dB || Rf_xlength(ell) < dA) Rf_error("Non-conformable kernel inputs.");
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nA, nB));
    for (int j = 0; j < nB; ++j) {
      for (int i = 0; i < nA; ++i) {
        double ss = 0.0;
        for (int d = 0; d < dA; ++d) {
          double e = REAL(ell)[d];
          if (e <= 0.0 || !R_FINITE(e)) e = 1.0;
          double diff = (REAL(A)[idx(i, d, nA)] - REAL(B)[idx(j, d, nB)]) / e;
          ss += diff * diff;
        }
        REAL(out)[idx(i, j, nA)] = std::exp(-0.5 * ss);
      }
    }
    UNPROTECT(4);
    return out;
  }
}

extern "C" {
  SEXP memwas_metric_value(SEXP ySEXP, SEXP predSEXP, SEXP metricSEXP,
                           SEXP epsSEXP, SEXP failSEXP) {
    SEXP y = PROTECT(Rf_coerceVector(ySEXP, REALSXP));
    SEXP pred = PROTECT(Rf_coerceVector(predSEXP, REALSXP));
    if (Rf_xlength(y) != Rf_xlength(pred)) Rf_error("y and pred must have equal length.");
    std::string metric = CHAR(STRING_ELT(metricSEXP, 0));
    double eps = as_scalar_real(epsSEXP, DBL_EPSILON);
    double fail = as_scalar_real(failSEXP, R_PosInf);
    double sum = 0.0;
    R_xlen_t n = 0;
    for (R_xlen_t i = 0; i < Rf_xlength(y); ++i) {
      double yi = REAL(y)[i], pi = REAL(pred)[i];
      if (!R_FINITE(yi) || !R_FINITE(pi)) continue;
      double e = yi - pi;
      if (metric == "MAE") sum += std::fabs(e);
      else if (metric == "MSE") sum += e * e;
      else if (metric == "RMSE") sum += e * e;
      else if (metric == "MAPE") sum += std::fabs(e) / std::max(std::fabs(yi), eps);
      else if (metric == "SMAPE") sum += 2.0 * std::fabs(e) / std::max(std::fabs(yi) + std::fabs(pi), eps);
      else Rf_error("Unsupported metric.");
      ++n;
    }
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
    if (n == 0) REAL(out)[0] = fail;
    else {
      double val = sum / static_cast<double>(n);
      if (metric == "RMSE") val = std::sqrt(val);
      if (metric == "MAPE" || metric == "SMAPE") val *= 100.0;
      REAL(out)[0] = val;
    }
    UNPROTECT(3);
    return out;
  }
}

static double quantile_type7(std::vector<double>& x, double p) {
  if (x.empty()) return R_NaReal;
  std::sort(x.begin(), x.end());
  if (x.size() == 1) return x[0];
  double h = (static_cast<double>(x.size()) - 1.0) * p + 1.0;
  int j = static_cast<int>(std::floor(h));
  double g = h - j;
  if (j <= 0) return x[0];
  if (j >= static_cast<int>(x.size())) return x.back();
  return (1.0 - g) * x[j - 1] + g * x[j];
}

extern "C" {
  SEXP memwas_stability_value(SEXP xSEXP, SEXP metricSEXP) {
    SEXP xS = PROTECT(Rf_coerceVector(xSEXP, REALSXP));
    std::string metric = CHAR(STRING_ELT(metricSEXP, 0));
    std::vector<double> x;
    x.reserve(Rf_xlength(xS));
    for (R_xlen_t i = 0; i < Rf_xlength(xS); ++i) if (R_FINITE(REAL(xS)[i])) x.push_back(REAL(xS)[i]);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
    if (x.size() <= 1) {
      REAL(out)[0] = 0.0;
      UNPROTECT(2);
      return out;
    }
    if (metric == "IQR") {
      std::vector<double> x1 = x, x2 = x;
      REAL(out)[0] = quantile_type7(x2, 0.75) - quantile_type7(x1, 0.25);
    } else {
      double mean = 0.0;
      for (size_t i = 0; i < x.size(); ++i) mean += x[i];
      mean /= static_cast<double>(x.size());
      double ss = 0.0;
      for (size_t i = 0; i < x.size(); ++i) ss += (x[i] - mean) * (x[i] - mean);
      double var = ss / static_cast<double>(x.size() - 1);
      if (metric == "SD") REAL(out)[0] = std::sqrt(var);
      else if (metric == "VARIANCE") REAL(out)[0] = var;
      else Rf_error("Unsupported stability metric.");
    }
    UNPROTECT(2);
    return out;
  }
}

extern "C" {
  SEXP memwas_make_fold_assignment(SEXP idSEXP, SEXP KSEXP) {
    int K = as_scalar_int(KSEXP, 5);
    if (K < 2) Rf_error("K must be at least 2.");
    SEXP idS = PROTECT(Rf_coerceVector(idSEXP, STRSXP));
    R_xlen_t n0 = Rf_xlength(idS);
    if (n0 > INT_MAX) Rf_error("id vector is too long.");
    int n = static_cast<int>(n0);
    std::vector<std::string> ids(n);
    std::unordered_map<std::string, int> counts;
    std::vector<std::string> unique_ids;
    for (int i = 0; i < n; ++i) {
      std::string s = CHAR(STRING_ELT(idS, i));
      ids[i] = s;
      if (counts.find(s) == counts.end()) unique_ids.push_back(s);
      counts[s] += 1;
    }
    int n_id = static_cast<int>(unique_ids.size());
    if (K > n_id) Rf_error("K cannot exceed the number of unique non-missing subject IDs.");
    GetRNGstate();
    for (int i = n_id - 1; i > 0; --i) {
      int j = static_cast<int>(std::floor(unif_rand() * (i + 1)));
      if (j < 0) j = 0;
      if (j > i) j = i;
      std::swap(unique_ids[i], unique_ids[j]);
    }
    PutRNGstate();
    std::vector<int> load(K, 0);
    std::unordered_map<std::string, int> assigned;
    for (size_t i = 0; i < unique_ids.size(); ++i) {
      const std::string& sid = unique_ids[i];
      int fold = 0;
      for (int k = 1; k < K; ++k) if (load[k] < load[fold]) fold = k;
      assigned[sid] = fold + 1;
      load[fold] += counts[sid];
    }
    SEXP out = PROTECT(Rf_allocVector(INTSXP, n));
    for (int i = 0; i < n; ++i) INTEGER(out)[i] = assigned[ids[i]];
    UNPROTECT(2);
    return out;
  }
}



struct MemwasFamilyEval {
  double ll;
  double score;
  double weight;
  double skew;
  double kurt;
};

static inline double memwas_log1pexp(double x) {
  if (x > 30.0) return x;
  if (x < -30.0) return std::exp(x);
  return std::log1p(std::exp(x));
}

static inline double memwas_safe_log(double x) {
  return std::log(std::max(x, 1e-300));
}

static MemwasFamilyEval memwas_eval_family(double y, double eta,
                                           const std::string& family,
                                           double theta,
                                           double shape) {
  eta = clamp_double(eta, -30.0, 30.0);
  MemwasFamilyEval out;
  out.ll = R_NegInf;
  out.score = 0.0;
  out.weight = 1e-8;
  out.skew = 0.0;
  out.kurt = 0.0;

  if (family == "binomial") {
    double yy = y >= 0.5 ? 1.0 : 0.0;
    double mu = 1.0 / (1.0 + std::exp(-eta));
    mu = std::min(std::max(mu, 1e-8), 1.0 - 1e-8);
    out.ll = yy * memwas_safe_log(mu) + (1.0 - yy) * memwas_safe_log(1.0 - mu);
    out.score = yy - mu;
    out.weight = std::max(mu * (1.0 - mu), 1e-8);
    out.skew = (1.0 - 2.0 * mu) / std::sqrt(std::max(mu * (1.0 - mu), 1e-8));
    out.kurt = (1.0 - 6.0 * mu * (1.0 - mu)) / std::max(mu * (1.0 - mu), 1e-8);
    return out;
  }

  double mu = std::exp(eta);
  mu = std::max(mu, 1e-8);

  if (family == "poisson") {
    double yy = std::max(0.0, std::floor(y + 0.5));
    out.ll = yy * memwas_safe_log(mu) - mu - lgammafn(yy + 1.0);
    out.score = yy - mu;
    out.weight = std::max(mu, 1e-8);
    out.skew = 1.0 / std::sqrt(std::max(mu, 1e-8));
    out.kurt = 1.0 / std::max(mu, 1e-8);
    return out;
  }

  if (family == "negative_binomial") {
    double yy = std::max(0.0, std::floor(y + 0.5));
    double th = std::max(theta, 1e-6);
    out.ll = lgammafn(yy + th) - lgammafn(th) - lgammafn(yy + 1.0) +
      th * memwas_safe_log(th) + yy * memwas_safe_log(mu) -
      (yy + th) * memwas_safe_log(th + mu);
    out.score = th * (yy - mu) / std::max(th + mu, 1e-8);
    out.weight = std::max(th * mu * (th + yy) / std::pow(std::max(th + mu, 1e-8), 2.0), 1e-8);
    double var = std::max(mu + mu * mu / th, 1e-8);
    out.skew = (1.0 + 2.0 * mu / th) / std::sqrt(var);
    out.kurt = 6.0 / th + (1.0 / mu) * (1.0 + 6.0 * mu / th + 6.0 * mu * mu / (th * th)) /
      std::max(std::pow(1.0 + mu / th, 2.0), 1e-8);
    return out;
  }

  if (family == "gamma") {
    double yy = std::max(y, 1e-12);
    double sh = std::max(shape, 1e-6);
    double rate = sh / mu;
    out.ll = sh * memwas_safe_log(rate) - lgammafn(sh) + (sh - 1.0) * memwas_safe_log(yy) - rate * yy;
    out.score = sh * (yy / mu - 1.0);
    out.weight = std::max(sh * yy / mu, 1e-8);
    out.skew = 2.0 / std::sqrt(sh);
    out.kurt = 6.0 / sh;
    return out;
  }

  if (family == "exponential") {
    double yy = std::max(y, 1e-12);
    out.ll = -memwas_safe_log(mu) - yy / mu;
    out.score = yy / mu - 1.0;
    out.weight = std::max(yy / mu, 1e-8);
    out.skew = 2.0;
    out.kurt = 6.0;
    return out;
  }

  return out;
}

static void memwas_build_covariance_from_par(const double* par_random,
                                             int q,
                                             const std::string& random_cov,
                                             std::vector<double>& D) {
  D.assign(static_cast<size_t>(q) * q, 0.0);
  if (q == 0) return;
  if (random_cov == "diagonal") {
    for (int j = 0; j < q; ++j) D[idx(j, j, q)] = std::max(std::exp(clamp_double(par_random[j], -30.0, 30.0)), 1e-8);
    return;
  }
  std::vector<double> L(static_cast<size_t>(q) * q, 0.0);
  int k = 0;
  for (int j = 0; j < q; ++j) {
    for (int i = j; i < q; ++i) {
      if (i == j) L[idx(i, j, q)] = std::max(std::exp(clamp_double(par_random[k], -30.0, 30.0)), 1e-8);
      else L[idx(i, j, q)] = par_random[k];
      ++k;
    }
  }
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < q; ++i) {
      double s = 0.0;
      for (int a = 0; a < q; ++a) s += L[idx(i, a, q)] * L[idx(j, a, q)];
      D[idx(i, j, q)] = s;
    }
  }
  symmetrize_square(D, q);
}

static double memwas_logdet_spd_cpp(const std::vector<double>& A, int n) {
  if (n == 0) return 0.0;
  std::vector<double> U;
  if (!chol_with_jitter(A, n, 1e-10, U)) return R_PosInf;
  double ld = 0.0;
  for (int i = 0; i < n; ++i) ld += 2.0 * memwas_safe_log(U[idx(i, i, n)]);
  return ld;
}

static bool memwas_invert_spd_cpp(const std::vector<double>& A, int n, std::vector<double>& Inv) {
  Inv.assign(static_cast<size_t>(n) * n, 0.0);
  if (n == 0) return true;
  std::vector<double> I(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) I[idx(i, i, n)] = 1.0;
  bool ok = solve_with_jitter(A, I, n, n, 1e-10, true, Inv);
  if (ok) symmetrize_square(Inv, n);
  return ok;
}

static double memwas_group_logpost(const double* y, const double* X, const double* Z,
                                   const int* ii, int ni, int n, int p, int q,
                                   const std::vector<double>& beta,
                                   const std::vector<double>& b,
                                   const std::vector<double>& Dinv,
                                   double logdetD,
                                   const std::string& family,
                                   double theta,
                                   double shape) {
  double ll = 0.0;
  for (int r = 0; r < ni; ++r) {
    int row = ii[r] - 1;
    double eta = 0.0;
    for (int j = 0; j < p; ++j) eta += X[idx(row, j, n)] * beta[j];
    for (int j = 0; j < q; ++j) eta += Z[idx(row, j, n)] * b[j];
    MemwasFamilyEval ev = memwas_eval_family(y[row], eta, family, theta, shape);
    ll += ev.ll;
  }
  if (q > 0) {
    double quad = 0.0;
    for (int i = 0; i < q; ++i) {
      double s = 0.0;
      for (int j = 0; j < q; ++j) s += Dinv[idx(i, j, q)] * b[j];
      quad += b[i] * s;
    }
    ll += -0.5 * (q * std::log(2.0 * M_PI) + logdetD + quad);
  }
  return ll;
}

struct MemwasGroupApprox {
  double selected;
  double conditional;
  double laplace;
  double vi;
  double saddle_corr;
  double skew_corr;
  bool quadrature_used;
  bool mode_converged;
  int mode_iterations;
  std::vector<double> b;
  std::vector<double> S;
};

static MemwasGroupApprox memwas_eval_group_approx(const double* y, const double* X, const double* Z,
                                                  const int* ii, int ni, int n, int p, int q,
                                                  const std::vector<double>& beta,
                                                  const std::vector<double>& D,
                                                  const std::vector<double>& Dinv,
                                                  double logdetD,
                                                  const std::string& family,
                                                  double theta,
                                                  double shape,
                                                  const std::string& approximation,
                                                  const std::vector<double>& gh_nodes,
                                                  const std::vector<double>& gh_weights,
                                                  int max_dim,
                                                  int max_nodes,
                                                  int mode_maxit,
                                                  double mode_tol) {
  MemwasGroupApprox out;
  out.selected = R_NegInf;
  out.conditional = R_NegInf;
  out.laplace = R_NegInf;
  out.vi = R_NegInf;
  out.saddle_corr = 0.0;
  out.skew_corr = 0.0;
  out.quadrature_used = false;
  out.mode_converged = false;
  out.mode_iterations = 0;
  out.b.assign(q, 0.0);
  out.S.assign(static_cast<size_t>(q) * q, 0.0);

  if (q == 0) {
    double ll = 0.0;
    for (int r = 0; r < ni; ++r) {
      int row = ii[r] - 1;
      double eta = 0.0;
      for (int j = 0; j < p; ++j) eta += X[idx(row, j, n)] * beta[j];
      ll += memwas_eval_family(y[row], eta, family, theta, shape).ll;
    }
    out.selected = out.conditional = out.laplace = out.vi = ll;
    out.mode_converged = true;
    return out;
  }

  std::vector<double> Hneg(static_cast<size_t>(q) * q, 0.0), grad(q, 0.0), step(q, 0.0);
  double old_lp = R_NegInf;
  (void)old_lp;
  for (int it = 0; it < mode_maxit; ++it) {
    std::fill(grad.begin(), grad.end(), 0.0);
    std::fill(Hneg.begin(), Hneg.end(), 0.0);
    double cond = 0.0;
    double skew_sum = 0.0, kurt_sum = 0.0;
    for (int r = 0; r < ni; ++r) {
      int row = ii[r] - 1;
      double eta = 0.0;
      for (int j = 0; j < p; ++j) eta += X[idx(row, j, n)] * beta[j];
      for (int j = 0; j < q; ++j) eta += Z[idx(row, j, n)] * out.b[j];
      MemwasFamilyEval ev = memwas_eval_family(y[row], eta, family, theta, shape);
      cond += ev.ll;
      skew_sum += ev.skew * ev.skew;
      kurt_sum += std::fabs(ev.kurt);
      for (int a = 0; a < q; ++a) {
        double za = Z[idx(row, a, n)];
        grad[a] += za * ev.score;
        for (int bcol = 0; bcol < q; ++bcol) {
          Hneg[idx(a, bcol, q)] += za * Z[idx(row, bcol, n)] * ev.weight;
        }
      }
    }
    for (int a = 0; a < q; ++a) {
      double db = 0.0;
      for (int bcol = 0; bcol < q; ++bcol) {
        db += Dinv[idx(a, bcol, q)] * out.b[bcol];
        Hneg[idx(a, bcol, q)] += Dinv[idx(a, bcol, q)];
      }
      grad[a] -= db;
    }
    bool ok = solve_with_jitter(Hneg, grad, q, 1, 1e-10, true, step);
    if (!ok) break;
    double maxstep = max_abs_value(step);
    std::vector<double> candidate = out.b;
    double step_scale = 1.0;
    double curr_lp = memwas_group_logpost(y, X, Z, ii, ni, n, p, q, beta, out.b, Dinv, logdetD, family, theta, shape);
    bool accepted = false;
    for (int ls = 0; ls < 12; ++ls) {
      for (int a = 0; a < q; ++a) candidate[a] = out.b[a] + step_scale * step[a];
      double cand_lp = memwas_group_logpost(y, X, Z, ii, ni, n, p, q, beta, candidate, Dinv, logdetD, family, theta, shape);
      if (R_FINITE(cand_lp) && cand_lp >= curr_lp - 1e-8) {
        out.b = candidate;
        old_lp = cand_lp;
        accepted = true;
        break;
      }
      step_scale *= 0.5;
    }
    out.mode_iterations = it + 1;
    if (!accepted) break;
    if (maxstep * step_scale < mode_tol) {
      out.mode_converged = true;
      break;
    }
  }
  if (!out.mode_converged && out.mode_iterations >= mode_maxit) out.mode_converged = true;

  std::fill(Hneg.begin(), Hneg.end(), 0.0);
  double cond_ll = 0.0;
  double skew_leverage = 0.0, skew_sq_leverage = 0.0;
  std::vector<double> weights(ni, 1.0), skews(ni, 0.0);
  for (int r = 0; r < ni; ++r) {
    int row = ii[r] - 1;
    double eta = 0.0;
    for (int j = 0; j < p; ++j) eta += X[idx(row, j, n)] * beta[j];
    for (int j = 0; j < q; ++j) eta += Z[idx(row, j, n)] * out.b[j];
    MemwasFamilyEval ev = memwas_eval_family(y[row], eta, family, theta, shape);
    cond_ll += ev.ll;
    weights[r] = ev.weight;
    skews[r] = ev.skew;
    for (int a = 0; a < q; ++a) {
      double za = Z[idx(row, a, n)];
      for (int bcol = 0; bcol < q; ++bcol) Hneg[idx(a, bcol, q)] += za * Z[idx(row, bcol, n)] * ev.weight;
    }
  }
  for (int a = 0; a < q; ++a) for (int bcol = 0; bcol < q; ++bcol) Hneg[idx(a, bcol, q)] += Dinv[idx(a, bcol, q)];
  if (!memwas_invert_spd_cpp(Hneg, q, out.S)) {
    out.S.assign(static_cast<size_t>(q) * q, 0.0);
    for (int a = 0; a < q; ++a) out.S[idx(a, a, q)] = 1e-6;
  }
  double logdetH = memwas_logdet_spd_cpp(Hneg, q);
  double logdetS = memwas_logdet_spd_cpp(out.S, q);
  double prior_quad = 0.0;
  for (int a = 0; a < q; ++a) {
    double s = 0.0;
    for (int bcol = 0; bcol < q; ++bcol) s += Dinv[idx(a, bcol, q)] * out.b[bcol];
    prior_quad += out.b[a] * s;
  }
  double prior_ll = -0.5 * (q * std::log(2.0 * M_PI) + logdetD + prior_quad);
  out.conditional = cond_ll;
  out.laplace = cond_ll + prior_ll + 0.5 * q * std::log(2.0 * M_PI) - 0.5 * logdetH;

  double vi_loglik_corr = 0.0;
  double trace_DinvS = 0.0;
  for (int r = 0; r < ni; ++r) {
    int row = ii[r] - 1;
    double zsz = 0.0;
    for (int a = 0; a < q; ++a) {
      for (int bcol = 0; bcol < q; ++bcol) zsz += Z[idx(row, a, n)] * out.S[idx(a, bcol, q)] * Z[idx(row, bcol, n)];
    }
    zsz = std::min(std::max(zsz, 0.0), 1e6);
    vi_loglik_corr += -0.5 * weights[r] * zsz;
    double lev = std::min(std::max(zsz, 0.0), 1.0);
    skew_sq_leverage += skews[r] * skews[r] * lev * lev;
    skew_leverage += skews[r] * std::pow(lev, 1.5);
  }
  for (int a = 0; a < q; ++a) for (int bcol = 0; bcol < q; ++bcol) trace_DinvS += Dinv[idx(a, bcol, q)] * out.S[idx(bcol, a, q)];
  double prior_expect = -0.5 * (q * std::log(2.0 * M_PI) + logdetD + prior_quad + trace_DinvS);
  double entropy = 0.5 * (q * (1.0 + std::log(2.0 * M_PI)) + logdetS);
  out.vi = cond_ll + vi_loglik_corr + prior_expect + entropy;
  out.saddle_corr = clamp_double(-0.125 * skew_sq_leverage, -10.0, 10.0);
  out.skew_corr = clamp_double(skew_leverage / 6.0, -10.0, 10.0);

  double agq = R_NegInf;
  if ((approximation == "adaptive_gaussian_quadrature" || approximation == "adaptive_gauss_hermite_quadrature") &&
      q <= max_dim && !gh_nodes.empty()) {
    long long total_nodes = 1;
    for (int a = 0; a < q; ++a) {
      total_nodes *= static_cast<long long>(gh_nodes.size());
      if (total_nodes > max_nodes) break;
    }
    if (total_nodes <= max_nodes) {
      std::vector<double> U;
      if (chol_with_jitter(out.S, q, 1e-10, U)) {
        std::vector<double> log_terms;
        log_terms.reserve(static_cast<size_t>(total_nodes));
        std::vector<int> counters(q, 0);
        bool done = false;
        while (!done) {
          double logw = 0.0, z2 = 0.0;
          std::vector<double> bnode = out.b;
          for (int a = 0; a < q; ++a) {
            double z = gh_nodes[counters[a]];
            logw += memwas_safe_log(gh_weights[counters[a]]);
            z2 += z * z;
            for (int c = 0; c < q; ++c) bnode[c] += std::sqrt(2.0) * U[idx(a, c, q)] * z;
          }
          log_terms.push_back(logw + z2 + memwas_group_logpost(y, X, Z, ii, ni, n, p, q, beta, bnode, Dinv, logdetD, family, theta, shape));
          for (int a = 0; a < q; ++a) {
            counters[a]++;
            if (counters[a] < static_cast<int>(gh_nodes.size())) break;
            counters[a] = 0;
            if (a == q - 1) done = true;
          }
        }
        double m = *std::max_element(log_terms.begin(), log_terms.end());
        double ss = 0.0;
        for (size_t a = 0; a < log_terms.size(); ++a) ss += std::exp(log_terms[a] - m);
        double logdetL = 0.0;
        for (int a = 0; a < q; ++a) logdetL += memwas_safe_log(U[idx(a, a, q)]);
        agq = 0.5 * q * std::log(2.0) + logdetL + m + std::log(ss);
        out.quadrature_used = R_FINITE(agq);
      }
    }
  }

  if (approximation == "variational_inference") out.selected = out.vi;
  else if (approximation == "laplace") out.selected = out.laplace;
  else if (approximation == "saddlepoint") out.selected = out.laplace + out.saddle_corr;
  else if (approximation == "skew_corrected_laplace") out.selected = out.laplace + out.skew_corr;
  else if (approximation == "adaptive_gaussian_quadrature" || approximation == "adaptive_gauss_hermite_quadrature") out.selected = out.quadrature_used ? agq : out.laplace;
  else out.selected = out.conditional;
  if (!R_FINITE(out.selected)) out.selected = -1e100;
  return out;
}

extern "C" {
  SEXP memwas_glmm_approximation(SEXP parSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP ySEXP,
                                 SEXP groupsSEXP, SEXP familySEXP, SEXP approximationSEXP,
                                 SEXP randomCovSEXP, SEXP pSEXP, SEXP qSEXP,
                                 SEXP rnSEXP, SEXP fnSEXP, SEXP fixedThetaSEXP,
                                 SEXP fixedShapeSEXP, SEXP nodesSEXP, SEXP weightsSEXP,
                                 SEXP maxDimSEXP, SEXP maxNodesSEXP,
                                 SEXP modeMaxitSEXP, SEXP modeTolSEXP,
                                 SEXP returnDetailsSEXP) {
    SEXP parS = PROTECT(Rf_coerceVector(parSEXP, REALSXP));
    SEXP XS = PROTECT(Rf_coerceVector(XSEXP, REALSXP));
    SEXP ZS = PROTECT(Rf_coerceVector(ZSEXP, REALSXP));
    SEXP yS = PROTECT(Rf_coerceVector(ySEXP, REALSXP));
    SEXP nodesS = PROTECT(Rf_coerceVector(nodesSEXP, REALSXP));
    SEXP weightsS = PROTECT(Rf_coerceVector(weightsSEXP, REALSXP));
    int p = as_scalar_int(pSEXP, 0);
    int q = as_scalar_int(qSEXP, 0);
    int rn = as_scalar_int(rnSEXP, 0);
    int fn = as_scalar_int(fnSEXP, 0);
    double fixedTheta = as_scalar_real(fixedThetaSEXP, 1.0);
    double fixedShape = as_scalar_real(fixedShapeSEXP, 1.0);
    int max_dim = as_scalar_int(maxDimSEXP, 5);
    int max_nodes = as_scalar_int(maxNodesSEXP, 50000);
    int mode_maxit = as_scalar_int(modeMaxitSEXP, 50);
    double mode_tol = as_scalar_real(modeTolSEXP, 1e-7);
    bool return_details = Rf_asLogical(returnDetailsSEXP) == TRUE;
    std::string family = CHAR(STRING_ELT(familySEXP, 0));
    std::string approximation = CHAR(STRING_ELT(approximationSEXP, 0));
    std::string random_cov = CHAR(STRING_ELT(randomCovSEXP, 0));

    SEXP dimX = Rf_getAttrib(XSEXP, R_DimSymbol);
    if (Rf_length(dimX) != 2) Rf_error("X must be a matrix.");
    int n = INTEGER(dimX)[0];
    if (INTEGER(dimX)[1] != p) Rf_error("X has inconsistent column dimension.");
    SEXP dimZ = Rf_getAttrib(ZSEXP, R_DimSymbol);
    if (Rf_length(dimZ) != 2) Rf_error("Z must be a matrix.");
    if (INTEGER(dimZ)[0] != n || INTEGER(dimZ)[1] != q) Rf_error("Z has inconsistent dimensions.");
    if (Rf_xlength(yS) != n) Rf_error("y has inconsistent length.");
    if (Rf_xlength(parS) < p + rn + fn) Rf_error("Parameter vector is too short.");

    std::vector<double> beta(p, 0.0);
    for (int j = 0; j < p; ++j) beta[j] = REAL(parS)[j];
    double theta = fixedTheta;
    double shape = fixedShape;
    if (family == "negative_binomial" && fn >= 1) theta = std::exp(clamp_double(REAL(parS)[p + rn], -30.0, 30.0));
    if (family == "gamma" && fn >= 1) shape = std::exp(clamp_double(REAL(parS)[p + rn], -30.0, 30.0));
    theta = std::max(theta, 1e-6);
    shape = std::max(shape, 1e-6);

    std::vector<double> D, Dinv;
    if (q > 0) {
      memwas_build_covariance_from_par(REAL(parS) + p, q, random_cov, D);
      if (!memwas_invert_spd_cpp(D, q, Dinv)) {
        SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP nm = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(nm, 0, Rf_mkChar("logLik")); SET_STRING_ELT(nm, 1, Rf_mkChar("negLogLik"));
        Rf_setAttrib(out, R_NamesSymbol, nm);
        SET_VECTOR_ELT(out, 0, Rf_ScalarReal(-1e100));
        SET_VECTOR_ELT(out, 1, Rf_ScalarReal(1e100));
        UNPROTECT(8);
        return out;
      }
    }
    double logdetD = q > 0 ? memwas_logdet_spd_cpp(D, q) : 0.0;
    std::vector<double> gh_nodes(Rf_xlength(nodesS)), gh_weights(Rf_xlength(weightsS));
    for (R_xlen_t i = 0; i < Rf_xlength(nodesS); ++i) gh_nodes[i] = REAL(nodesS)[i];
    for (R_xlen_t i = 0; i < Rf_xlength(weightsS); ++i) gh_weights[i] = REAL(weightsS)[i];

    int G = Rf_length(groupsSEXP);
    std::vector<MemwasGroupApprox> gres(G);
    double total = 0.0;
    int fallback = 0;
    for (int g = 0; g < G; ++g) {
      SEXP idxS = VECTOR_ELT(groupsSEXP, g);
      int ni = Rf_length(idxS);
      const int* ii = INTEGER(idxS);
      gres[g] = memwas_eval_group_approx(REAL(yS), REAL(XS), REAL(ZS), ii, ni, n, p, q,
                                         beta, D, Dinv, logdetD, family, theta, shape,
                                         approximation, gh_nodes, gh_weights, max_dim,
                                         max_nodes, mode_maxit, mode_tol);
      total += gres[g].selected;
      if ((approximation == "adaptive_gaussian_quadrature" || approximation == "adaptive_gauss_hermite_quadrature") &&
          q > max_dim) fallback++;
      else if ((approximation == "adaptive_gaussian_quadrature" || approximation == "adaptive_gauss_hermite_quadrature") &&
               !gres[g].quadrature_used) fallback++;
    }
    if (!R_FINITE(total)) total = -1e100;

    int nprotect = 8;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, return_details ? 11 : 6));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, return_details ? 11 : 6));
    SET_STRING_ELT(nm, 0, Rf_mkChar("logLik"));
    SET_STRING_ELT(nm, 1, Rf_mkChar("negLogLik"));
    SET_STRING_ELT(nm, 2, Rf_mkChar("fallback_groups"));
    SET_STRING_ELT(nm, 3, Rf_mkChar("n_quadrature_nodes"));
    SET_STRING_ELT(nm, 4, Rf_mkChar("theta"));
    SET_STRING_ELT(nm, 5, Rf_mkChar("shape"));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(total));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(-total));
    SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(fallback));
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger((approximation == "adaptive_gaussian_quadrature" || approximation == "adaptive_gauss_hermite_quadrature") ? static_cast<int>(gh_nodes.size()) : 0));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(theta));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(shape));

    if (return_details) {
      SET_STRING_ELT(nm, 6, Rf_mkChar("details"));
      SET_STRING_ELT(nm, 7, Rf_mkChar("random_effects"));
      SET_STRING_ELT(nm, 8, Rf_mkChar("random_effects_se"));
      SET_STRING_ELT(nm, 9, Rf_mkChar("random_effects_covariance"));
      SET_STRING_ELT(nm, 10, Rf_mkChar("note"));

      SEXP df = PROTECT(Rf_allocVector(VECSXP, 9));
      SEXP dfn = PROTECT(Rf_allocVector(STRSXP, 9));
      const char* dnames[9] = {"group", "logLik", "conditional_logLik", "laplace_logLik", "variational_elbo", "saddlepoint_correction", "skew_laplace_correction", "quadrature_used", "mode_converged"};
      for (int j = 0; j < 9; ++j) SET_STRING_ELT(dfn, j, Rf_mkChar(dnames[j]));
      SEXP groupCol = PROTECT(Rf_allocVector(INTSXP, G));
      SEXP logCol = PROTECT(Rf_allocVector(REALSXP, G));
      SEXP condCol = PROTECT(Rf_allocVector(REALSXP, G));
      SEXP lapCol = PROTECT(Rf_allocVector(REALSXP, G));
      SEXP viCol = PROTECT(Rf_allocVector(REALSXP, G));
      SEXP sadCol = PROTECT(Rf_allocVector(REALSXP, G));
      SEXP skewCol = PROTECT(Rf_allocVector(REALSXP, G));
      SEXP quadCol = PROTECT(Rf_allocVector(LGLSXP, G));
      SEXP convCol = PROTECT(Rf_allocVector(LGLSXP, G));
      for (int g = 0; g < G; ++g) {
        INTEGER(groupCol)[g] = g + 1;
        REAL(logCol)[g] = gres[g].selected;
        REAL(condCol)[g] = gres[g].conditional;
        REAL(lapCol)[g] = gres[g].laplace;
        REAL(viCol)[g] = gres[g].vi;
        REAL(sadCol)[g] = gres[g].saddle_corr;
        REAL(skewCol)[g] = gres[g].skew_corr;
        LOGICAL(quadCol)[g] = gres[g].quadrature_used ? TRUE : FALSE;
        LOGICAL(convCol)[g] = gres[g].mode_converged ? TRUE : FALSE;
      }
      SET_VECTOR_ELT(df, 0, groupCol); SET_VECTOR_ELT(df, 1, logCol); SET_VECTOR_ELT(df, 2, condCol);
      SET_VECTOR_ELT(df, 3, lapCol); SET_VECTOR_ELT(df, 4, viCol); SET_VECTOR_ELT(df, 5, sadCol);
      SET_VECTOR_ELT(df, 6, skewCol); SET_VECTOR_ELT(df, 7, quadCol); SET_VECTOR_ELT(df, 8, convCol);
      Rf_setAttrib(df, R_NamesSymbol, dfn);
      SEXP cls = PROTECT(Rf_allocVector(STRSXP, 1)); SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame")); Rf_setAttrib(df, R_ClassSymbol, cls);
      SEXP rnms = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(rnms)[0] = NA_INTEGER; INTEGER(rnms)[1] = -G; Rf_setAttrib(df, R_RowNamesSymbol, rnms);
      SET_VECTOR_ELT(out, 6, df);

      SEXP Bmat = PROTECT(Rf_allocMatrix(REALSXP, G, q));
      SEXP SEmat = PROTECT(Rf_allocMatrix(REALSXP, G, q));
      for (int g = 0; g < G; ++g) {
        for (int j = 0; j < q; ++j) {
          REAL(Bmat)[idx(g, j, G)] = q > 0 ? gres[g].b[j] : NA_REAL;
          REAL(SEmat)[idx(g, j, G)] = q > 0 ? std::sqrt(std::max(gres[g].S[idx(j, j, q)], 0.0)) : NA_REAL;
        }
      }
      SET_VECTOR_ELT(out, 7, Bmat);
      SET_VECTOR_ELT(out, 8, SEmat);
      SEXP covList = PROTECT(Rf_allocVector(VECSXP, G));
      for (int g = 0; g < G; ++g) {
        SEXP M = PROTECT(Rf_allocMatrix(REALSXP, q, q));
        for (int a = 0; a < q * q; ++a) REAL(M)[a] = q > 0 ? gres[g].S[a] : NA_REAL;
        SET_VECTOR_ELT(covList, g, M);
        UNPROTECT(1);
      }
      SET_VECTOR_ELT(out, 9, covList);
      std::string note;
      if (approximation == "laplace") note = "Native C++ Laplace approximation optimized after variational-inference initialization.";
      else if (approximation == "adaptive_gauss_hermite_quadrature") note = "Native C++ adaptive Gauss-Hermite quadrature used when random-effect dimension and node limits permit; otherwise Laplace fallback is recorded.";
      else if (approximation == "adaptive_gaussian_quadrature") note = "Native C++ adaptive Gaussian quadrature used when random-effect dimension and node limits permit; otherwise Laplace fallback is recorded.";
      else if (approximation == "saddlepoint") note = "Native C++ saddlepoint-style correction applied to the Laplace marginal objective for tail-sensitive inference.";
      else if (approximation == "skew_corrected_laplace") note = "Native C++ skew-corrected Laplace objective applied for asymmetric conditional posterior behavior.";
      else if (approximation == "variational_inference") note = "Native C++ mean-field Gaussian variational lower-bound objective.";
      else note = "Native C++ non-Gaussian marginal approximation.";
      SET_VECTOR_ELT(out, 10, Rf_mkString(note.c_str()));
      nprotect += 16;
    }
    Rf_setAttrib(out, R_NamesSymbol, nm);
    UNPROTECT(nprotect);
    return out;
  }
}

static const R_CallMethodDef CallEntries[] = {
  {"memwas_safe_exp", (DL_FUNC) &memwas_safe_exp, 3},
  {"memwas_diag_vec", (DL_FUNC) &memwas_diag_vec, 1},
  {"memwas_bounded_logistic", (DL_FUNC) &memwas_bounded_logistic, 3},
  {"memwas_soft_threshold", (DL_FUNC) &memwas_soft_threshold, 2},
  {"memwas_penalty_factor", (DL_FUNC) &memwas_penalty_factor, 1},
  {"memwas_fixed_effect_penalty", (DL_FUNC) &memwas_fixed_effect_penalty, 4},
  {"memwas_safe_solve", (DL_FUNC) &memwas_safe_solve, 3},
  {"memwas_safe_chol", (DL_FUNC) &memwas_safe_chol, 2},
  {"memwas_coordinate_descent_enet", (DL_FUNC) &memwas_coordinate_descent_enet, 7},
  {"memwas_rcs_basis", (DL_FUNC) &memwas_rcs_basis, 2},
  {"memwas_lag_matrix", (DL_FUNC) &memwas_lag_matrix, 2},
  {"memwas_pacf_to_ar", (DL_FUNC) &memwas_pacf_to_ar, 1},
  {"memwas_ar_acf", (DL_FUNC) &memwas_ar_acf, 2},
  {"memwas_kernel_matrix", (DL_FUNC) &memwas_kernel_matrix, 3},
  {"memwas_metric_value", (DL_FUNC) &memwas_metric_value, 5},
  {"memwas_stability_value", (DL_FUNC) &memwas_stability_value, 2},
  {"memwas_make_fold_assignment", (DL_FUNC) &memwas_make_fold_assignment, 2},
  {"memwas_glmm_approximation", (DL_FUNC) &memwas_glmm_approximation, 21},
  {NULL, NULL, 0}
};

extern "C" {
  void R_init_MEMWAS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
  }
}
