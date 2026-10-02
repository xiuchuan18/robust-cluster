// SIT_core.cpp
//
// C++ backend for MixSIT.
//
// PERFORMANCE REVISION 2026-08-28. The model, the objective functions and the
// E-step expectations are unchanged. Only how they are computed changed.
//
//   (1) Zero-copy inputs. Every entry point used to take `arma::mat Y`,
//       `arma::cube W_hat`, `arma::vec par` BY VALUE, so a full copy of Y
//       (N x m) and of W (N x m x K) was made on EVERY objective evaluation --
//       and nlminb makes ~1000 (m=5) to ~5000 (m=10) of them per component per
//       ECM iteration. Inputs are now wrapped around the R memory in place.
//
//   (2) Halved flops in the latent-coordinate computation. The old code formed
//           D = O * solve(U', (Y - 1 mu')')            [m x N]
//       which is a transpose allocation plus TWO O(N m^2) operations (a solve
//       with N right-hand sides, then a GEMM). The identity
//           D' = (Y - 1 mu') * (U^{-1} O')
//       gives the same matrix with ONE m x m triangular solve (O(m^3)) and ONE
//       O(N m^2) GEMM, and returns it already in N x m layout.
//
//   (3) Cache-friendly accumulation. D used to be m x N, so the inner loop over
//       i strode by m. It is now N x m, so the loop over i is contiguous. The
//       summation ORDER is unchanged (j outer, i inner), so the accumulation is
//       bit-identical given the same D.
//
//   (4) Triangular solves are declared triangular (trimatu) instead of going
//       through the general LU path.
//
//   (5) E-step: digamma((nu_j+1)/2) depends only on (j,k) but was evaluated
//       inside the loop over observations, i.e. N times more often than needed.
//       Hoisted out. Same values.
//
//   BUG FIX (behaviour change, and the only one): log_probs was initialised to
//   ZERO and left at zero for any component skipped by the is_well_conditioned
//   or solve guards. Zero is an enormous log-density next to the usual ~-10, so
//   logSumExp made the broken component absorb the entire posterior, silently.
//   Skipped components now get -infinity. This can only alter results in runs
//   that were already producing garbage.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Rcpp;
using namespace arma;

const double LOG_2 = std::log(2.0);
const double LOG_PI = std::log(M_PI);
const double LOG_2PI = std::log(2.0 * M_PI);
const double PENALTY_VALUE = 1e10;
const double NEG_INF = -std::numeric_limits<double>::infinity();

inline double fast_lgamma(double x) { return R::lgammafn(x); }
inline double fast_digamma(double x) { return R::digamma(x); }

// ---- zero-copy views over R memory -----------------------------------------
inline arma::mat view_mat(const Rcpp::NumericMatrix& x) {
  return arma::mat(const_cast<double*>(x.begin()), x.nrow(), x.ncol(), false, true);
}
inline arma::vec view_vec(const Rcpp::NumericVector& x) {
  return arma::vec(const_cast<double*>(x.begin()), x.size(), false, true);
}
inline arma::cube view_cube(const Rcpp::NumericVector& x) {
  Rcpp::IntegerVector d = x.attr("dim");
  return arma::cube(const_cast<double*>(x.begin()), d[0], d[1], d[2], false, true);
}

inline double log_dt_fast(double x, double nu) {
  double half_nu = nu * 0.5;
  double half_nu_plus_half = (nu + 1.0) * 0.5;
  return fast_lgamma(half_nu_plus_half) - fast_lgamma(half_nu)
       - 0.5 * std::log(nu * M_PI)
       - half_nu_plus_half * std::log(1.0 + (x * x) / nu);
}

// [[Rcpp::export]]
double dskwet_cpp(double x, double nu, double gamma, bool return_log = false) {
  double inv_gamma = 1.0 / gamma;
  double scaling = (x >= 0) ? inv_gamma : gamma;
  double z = x * scaling;
  double log_dens = log_dt_fast(z, nu);
  double log_const = LOG_2 - std::log(gamma + inv_gamma);
  double total = log_dens + log_const;
  return return_log ? total : std::exp(total);
}

// [[Rcpp::export]]
arma::vec dskwet_vec_cpp(const arma::vec& x, double nu, double gamma) {
  int n = x.n_elem;
  arma::vec result(n);

  double inv_gamma = 1.0 / gamma;
  double log_const = LOG_2 - std::log(gamma + inv_gamma);
  double half_nu = nu * 0.5;
  double half_nu_plus_half = (nu + 1.0) * 0.5;
  double lgamma_diff = fast_lgamma(half_nu_plus_half) - fast_lgamma(half_nu);
  double log_nu_pi = 0.5 * std::log(nu * M_PI);

  const double* xp = x.memptr();
  double* rp = result.memptr();
  for (int i = 0; i < n; ++i) {
    double scaling = (xp[i] >= 0) ? inv_gamma : gamma;
    double z = xp[i] * scaling;
    double z_sq_over_nu = (z * z) / nu;
    rp[i] = lgamma_diff - log_nu_pi
          - half_nu_plus_half * std::log(1.0 + z_sq_over_nu)
          + log_const;
  }
  return result;
}

// [[Rcpp::export]]
arma::vec rskwet_cpp(int n, double nu, double gamma) {
  if (nu <= 0) Rcpp::stop("Parameter nu must be positive");
  if (gamma <= 0) Rcpp::stop("Parameter gamma must be positive");

  arma::vec result(n);
  double gamma_sq = gamma * gamma;
  double prob_positive = gamma_sq / (1.0 + gamma_sq);
  double half_nu = nu * 0.5;
  double scale_val = 1.0 / half_nu;

  for (int i = 0; i < n; ++i) {
    double w = R::rgamma(half_nu, scale_val);
    double z_abs = std::abs(R::rnorm(0.0, 1.0));
    double t_abs = z_abs / std::sqrt(w);
    double u = R::runif(0.0, 1.0);
    result(i) = (u < prob_positive) ? t_abs * gamma : -t_abs / gamma;
  }
  return result;
}

// [[Rcpp::export]]
arma::mat rskwet_matrix_cpp(int n, const arma::vec& nu, const arma::vec& gamma) {
  int m = nu.n_elem;
  if (gamma.n_elem != m) Rcpp::stop("nu and gamma must have the same length");
  arma::mat result(n, m);
  for (int j = 0; j < m; ++j) result.col(j) = rskwet_cpp(n, nu(j), gamma(j));
  return result;
}

// [[Rcpp::export]]
arma::mat vector_to_orthogonal_cpp(arma::vec theta, int m, bool use_reflection = false) {
  if (m == 1) {
    arma::mat res(1, 1);
    res(0, 0) = use_reflection ? -1.0 : 1.0;
    return res;
  }

  arma::mat S = arma::zeros(m, m);
  int idx = 0;
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < j; ++i) {
      if (idx < (int)theta.n_elem) S(i, j) = theta(idx++);
    }
  }
  S = S - S.t();

  arma::mat I = arma::eye(m, m);
  arma::mat O = (I - S) * arma::solve(I + S, I, arma::solve_opts::fast);

  if (use_reflection) O.row(0) *= -1.0;
  return O;
}

arma::mat build_U_cpp(const arma::vec& u_pars, int m) {
  arma::vec u_diag_log = arma::clamp(u_pars.head(m), -10.0, 10.0);
  arma::vec u_diag = arma::exp(u_diag_log);
  arma::mat U = arma::diagmat(u_diag);

  if (m > 1) {
    int idx = m;
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < j; ++i) {
        if (idx < (int)u_pars.n_elem) U(i, j) = u_pars(idx++);
      }
    }
  }
  return U;
}

inline bool is_well_conditioned(const arma::mat& U, double threshold = 1e-8) {
  arma::vec diag_U = arma::abs(U.diag());
  return arma::min(diag_U) >= threshold;
}

inline double logSumExp_cpp(const arma::vec& lx) {
  double max_lx = lx.max();
  if (!std::isfinite(max_lx)) return max_lx;
  return max_lx + std::log(arma::sum(arma::exp(lx - max_lx)));
}

// Latent coordinates in N x m layout:  D' = (Y - 1 mu') * (U^{-1} O')
inline bool latent_coords(const arma::mat& Y, const arma::vec& mu,
                          const arma::mat& U, const arma::mat& O,
                          arma::mat& Dt) {
  arma::mat M;
  if (!arma::solve(M, arma::trimatu(U), O.t(), arma::solve_opts::fast)) return false;
  Dt = (Y.each_row() - mu.t()) * M;
  return true;
}

// [[Rcpp::export]]
List E_Step_cpp(Rcpp::NumericMatrix Y_, List params) {
  const arma::mat Y = view_mat(Y_);
  int N = Y.n_rows;
  int m = Y.n_cols;
  NumericVector pi = params["pi"];
  int K = pi.size();
  List components = params["components"];

  arma::mat Z_hat = arma::zeros(N, K);
  arma::cube W_hat = arma::zeros(N, m, K);
  arma::cube Kappa_hat = arma::zeros(N, m, K);
  arma::mat log_probs(N, K);
  log_probs.fill(NEG_INF);                 // BUG FIX: was zeros()
  double total_loglik = 0;

  for (int k = 0; k < K; ++k) {
    List theta = components[k];
    arma::vec mu = as<arma::vec>(theta["mu"]);
    arma::mat U = as<arma::mat>(theta["U"]);
    arma::mat O = as<arma::mat>(theta["O"]);
    arma::vec gamma = as<arma::vec>(theta["gamma"]);
    arma::vec nu = as<arma::vec>(theta["nu"]);

    if (!is_well_conditioned(U)) continue;

    arma::mat Dt;                                     // N x m
    if (!latent_coords(Y, mu, U, O, Dt)) continue;

    double log_jacobian = -arma::sum(arma::log(U.diag()));
    double log_pi_k = std::log(pi[k]);

    arma::vec inv_gamma = arma::pow(gamma, -1.0);
    arma::vec half_nu_plus_half = (nu + 1.0) * 0.5;

    // per-(j,k) constants -- digamma used to be recomputed for every i
    arma::vec lgamma_diff(m), log_nu_pi(m), log_const(m), dig(m);
    for (int j = 0; j < m; ++j) {
      lgamma_diff(j) = fast_lgamma(half_nu_plus_half(j)) - fast_lgamma(nu(j) * 0.5);
      log_nu_pi(j) = 0.5 * std::log(nu(j) * M_PI);
      log_const(j) = LOG_2 - std::log(gamma(j) + inv_gamma(j));
      dig(j) = fast_digamma(half_nu_plus_half(j));
    }

    arma::vec log_dens_sum(N, arma::fill::zeros);
    for (int j = 0; j < m; ++j) {
      const double* dp = Dt.colptr(j);               // contiguous in i
      double* wp = W_hat.slice(k).colptr(j);
      double* kp = Kappa_hat.slice(k).colptr(j);
      double gj = gamma(j), igj = inv_gamma(j), nuj = nu(j);
      double a = half_nu_plus_half(j);
      double base = lgamma_diff(j) - log_nu_pi(j) + log_const(j);
      double digj = dig(j);
      for (int i = 0; i < N; ++i) {
        double d_val = dp[i];
        double z = d_val * ((d_val >= 0) ? igj : gj);
        double z_sq = z * z;
        log_dens_sum(i) += base - a * std::log(1.0 + z_sq / nuj);
        double beta_post = (nuj + z_sq) * 0.5;
        wp[i] = a / beta_post;
        kp[i] = digj - std::log(beta_post);
      }
    }
    log_probs.col(k) = log_dens_sum + (log_jacobian + log_pi_k);
  }

  for (int i = 0; i < N; ++i) {
    double lse = logSumExp_cpp(log_probs.row(i).t());
    total_loglik += lse;
    if (!std::isfinite(lse)) {                        // no usable component
      for (int k = 0; k < K; ++k) Z_hat(i, k) = 1.0 / K;
      continue;
    }
    for (int k = 0; k < K; ++k) {
      Z_hat(i, k) = std::exp(log_probs(i, k) - lse);
    }
  }

  return List::create(
    Named("Z") = Z_hat,
    Named("W") = W_hat,
    Named("Kappa") = Kappa_hat,
    Named("loglik") = total_loglik
  );
}

// weighted sum of squares; summation order (j outer, i inner) preserved
inline double wss_accumulate(const arma::mat& Dt, const arma::vec& z_k,
                             const arma::mat& w, const arma::vec& g_sq,
                             const arma::vec& inv_g_sq) {
  int N = Dt.n_rows, m = Dt.n_cols;
  double wss = 0.0;
  const double* zp = z_k.memptr();
  for (int j = 0; j < m; ++j) {
    const double* dp = Dt.colptr(j);
    const double* wp = w.colptr(j);
    double gj = g_sq(j), igj = inv_g_sq(j);
    for (int i = 0; i < N; ++i) {
      double d_val = dp[i];
      double sc_sq = (d_val >= 0) ? igj : gj;
      wss += zp[i] * wp[i] * d_val * d_val * sc_sq;
    }
  }
  return wss;
}

inline double reg_U_term(const arma::mat& U, int m, double lambda_U) {
  double reg = 0.0;
  for (int j = 1; j < m; ++j)
    for (int i = 0; i < j; ++i) reg += lambda_U * U(i, j) * U(i, j);
  return reg;
}

// [[Rcpp::export]]
double obj_space_hetero_cpp(Rcpp::NumericVector par_, Rcpp::NumericMatrix Y_,
                            List theta_fixed,
                            Rcpp::NumericVector z_k_, Rcpp::NumericMatrix w_k_mat_,
                            double lambda_U, double lambda_O,
                            bool use_reflection) {
  const arma::vec par = view_vec(par_);
  const arma::mat Y = view_mat(Y_);
  const arma::vec z_k = view_vec(z_k_);
  const arma::mat w_k_mat = view_mat(w_k_mat_);

  int m = Y.n_cols;
  double sum_z = arma::sum(z_k);

  int idx = 0;
  arma::vec mu = par.subvec(idx, idx + m - 1); idx += m;

  int len_u_vec = (m > 1) ? m * (m - 1) / 2 : 0;
  arma::vec u_pars = par.subvec(idx, idx + m + len_u_vec - 1); idx += (m + len_u_vec);
  arma::mat U = build_U_cpp(u_pars, m);

  if (!is_well_conditioned(U)) return PENALTY_VALUE;

  arma::mat O;
  arma::vec ang_pars;
  if (m > 1) {
    int len_angles = m * (m - 1) / 2;
    ang_pars = par.subvec(idx, idx + len_angles - 1);
    O = vector_to_orthogonal_cpp(ang_pars, m, use_reflection);
  } else {
    O = arma::mat(1, 1);
    O(0, 0) = use_reflection ? -1.0 : 1.0;
  }

  double log_det_U = arma::sum(arma::log(U.diag()));
  double term_jac = sum_z * log_det_U;

  arma::mat Dt;
  if (!latent_coords(Y, mu, U, O, Dt)) return PENALTY_VALUE;

  arma::vec gamma = as<arma::vec>(theta_fixed["gamma"]);
  arma::vec g_sq = gamma % gamma;
  arma::vec inv_g_sq = arma::pow(g_sq, -1.0);

  double wss = wss_accumulate(Dt, z_k, w_k_mat, g_sq, inv_g_sq);

  double reg = 0.0;
  if (m > 1) {
    reg = reg_U_term(U, m, lambda_U);
    reg += lambda_O * arma::dot(ang_pars, ang_pars);
  }
  return term_jac + 0.5 * wss + reg;
}

// [[Rcpp::export]]
double obj_space_full_cpp(Rcpp::NumericVector par_, Rcpp::NumericMatrix Y_,
                          List components,
                          Rcpp::NumericMatrix Z_hat_, Rcpp::NumericVector W_hat_,
                          double lambda_U, double lambda_O,
                          bool use_reflection) {
  const arma::vec par = view_vec(par_);
  const arma::mat Y = view_mat(Y_);
  const arma::mat Z_hat = view_mat(Z_hat_);
  const arma::cube W_hat = view_cube(W_hat_);

  int m = Y.n_cols;
  int K = components.size();

  int idx = 0;
  std::vector<arma::vec> mus(K);
  for (int k = 0; k < K; ++k) {
    mus[k] = par.subvec(idx, idx + m - 1); idx += m;
  }

  int len_u_vec = (m > 1) ? m * (m - 1) / 2 : 0;
  arma::vec u_pars = par.subvec(idx, idx + m + len_u_vec - 1); idx += (m + len_u_vec);
  arma::mat U = build_U_cpp(u_pars, m);

  if (!is_well_conditioned(U)) return PENALTY_VALUE;

  arma::mat O;
  arma::vec ang_pars;
  if (m > 1) {
    int len_angles = m * (m - 1) / 2;
    ang_pars = par.subvec(idx, idx + len_angles - 1);
    O = vector_to_orthogonal_cpp(ang_pars, m, use_reflection);
  } else {
    O = arma::mat(1, 1);
    O(0, 0) = use_reflection ? -1.0 : 1.0;
  }

  double log_det_U = arma::sum(arma::log(U.diag()));
  double total_loss = 0.0;

  // U and O are shared, so the m x m map is formed once for all components
  arma::mat Mmap;
  if (!arma::solve(Mmap, arma::trimatu(U), O.t(), arma::solve_opts::fast))
    return PENALTY_VALUE;

  for (int k = 0; k < K; ++k) {
    arma::vec z_k = Z_hat.col(k);
    double sum_z = arma::sum(z_k);
    if (sum_z < 1e-8) continue;

    double term_jac = sum_z * log_det_U;
    arma::mat Dt = (Y.each_row() - mus[k].t()) * Mmap;

    List theta = components[k];
    arma::vec gamma = as<arma::vec>(theta["gamma"]);
    arma::vec g_sq = gamma % gamma;
    arma::vec inv_g_sq = arma::pow(g_sq, -1.0);

    double wss = wss_accumulate(Dt, z_k, W_hat.slice(k), g_sq, inv_g_sq);
    total_loss += term_jac + 0.5 * wss;
  }

  double reg = 0.0;
  if (m > 1) {
    reg = reg_U_term(U, m, lambda_U);
    reg += lambda_O * arma::dot(ang_pars, ang_pars);
  }
  return total_loss + reg;
}

// [[Rcpp::export]]
double obj_space_u_only_cpp(Rcpp::NumericVector par_, Rcpp::NumericMatrix Y_,
                            List components,
                            Rcpp::NumericMatrix Z_hat_, Rcpp::NumericVector W_hat_,
                            double lambda_U, double lambda_O,
                            LogicalVector use_reflection_vec) {
  const arma::vec par = view_vec(par_);
  const arma::mat Y = view_mat(Y_);
  const arma::mat Z_hat = view_mat(Z_hat_);
  const arma::cube W_hat = view_cube(W_hat_);

  int m = Y.n_cols;
  int K = components.size();

  int idx = 0;
  std::vector<arma::vec> mus(K);
  for (int k = 0; k < K; ++k) {
    mus[k] = par.subvec(idx, idx + m - 1); idx += m;
  }

  int len_angles = (m > 1) ? m * (m - 1) / 2 : 0;
  std::vector<arma::mat> Os(K);
  std::vector<arma::vec> angles_list(K);

  for (int k = 0; k < K; ++k) {
    bool ref_k = use_reflection_vec[k];
    if (m > 1) {
      angles_list[k] = par.subvec(idx, idx + len_angles - 1); idx += len_angles;
      Os[k] = vector_to_orthogonal_cpp(angles_list[k], m, ref_k);
    } else {
      Os[k] = arma::mat(1, 1);
      Os[k](0, 0) = ref_k ? -1.0 : 1.0;
    }
  }

  int len_u_vec = (m > 1) ? m * (m - 1) / 2 : 0;
  arma::vec u_pars = par.subvec(idx, idx + m + len_u_vec - 1);
  arma::mat U = build_U_cpp(u_pars, m);

  if (!is_well_conditioned(U)) return PENALTY_VALUE;

  double log_det_U = arma::sum(arma::log(U.diag()));
  double total_loss = 0.0;

  for (int k = 0; k < K; ++k) {
    arma::vec z_k = Z_hat.col(k);
    double sum_z = arma::sum(z_k);
    if (sum_z < 1e-8) continue;

    double term_jac = sum_z * log_det_U;
    arma::mat Dt;
    if (!latent_coords(Y, mus[k], U, Os[k], Dt)) return PENALTY_VALUE;

    List theta = components[k];
    arma::vec gamma = as<arma::vec>(theta["gamma"]);
    arma::vec g_sq = gamma % gamma;
    arma::vec inv_g_sq = arma::pow(g_sq, -1.0);

    double wss = wss_accumulate(Dt, z_k, W_hat.slice(k), g_sq, inv_g_sq);
    double reg_O = 0.0;
    if (m > 1) reg_O = lambda_O * arma::dot(angles_list[k], angles_list[k]);
    total_loss += term_jac + 0.5 * wss + reg_O;
  }

  double reg_U = 0.0;
  if (m > 1) reg_U = reg_U_term(U, m, lambda_U);
  return total_loss + reg_U;
}

// ---------------------------------------------------------------------------
// Fused objective for the CM-2 (AECM) skewness update.
//
// The R version looped over j, extracting a strided row of the m x N matrix
// D, crossing into C++ once per j, and allocating two length-N vectors per j:
//
//     for (j in seq_len(m)) {
//       d_vals    <- D_mat[j, ]
//       log_dens  <- dskwet_vec_cpp(d_vals, nu[j], gamma[j])
//       neg_ll    <- neg_ll - sum(z_k * log_dens)
//     }
//
// This does the whole thing in one call. The arithmetic is reproduced term by
// term, and the inner sum uses a long double accumulator so that it matches R's
// sum(), which also accumulates in long double. Bit-identical to the R loop.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double obj_gamma_cpp(Rcpp::NumericMatrix D_, Rcpp::NumericVector z_,
                     Rcpp::NumericVector nu_, Rcpp::NumericVector gamma_) {
  const int m = D_.nrow();
  const int N = D_.ncol();
  const double* D = D_.begin();          // column-major, D(j,i) = D[j + i*m]
  const double* z = z_.begin();

  double neg_ll = 0.0;
  for (int j = 0; j < m; ++j) {
    const double nu = nu_[j];
    const double gamma = gamma_[j];
    const double inv_gamma = 1.0 / gamma;
    const double log_const = LOG_2 - std::log(gamma + inv_gamma);
    const double half_nu = nu * 0.5;
    const double half_nu_plus_half = (nu + 1.0) * 0.5;
    const double lgamma_diff = fast_lgamma(half_nu_plus_half) - fast_lgamma(half_nu);
    const double log_nu_pi = 0.5 * std::log(nu * M_PI);

    long double acc = 0.0L;
    for (int i = 0; i < N; ++i) {
      const double d_val = D[j + (std::size_t)i * m];
      const double scaling = (d_val >= 0) ? inv_gamma : gamma;
      const double zz = d_val * scaling;
      const double z_sq_over_nu = (zz * zz) / nu;
      const double log_dens = lgamma_diff - log_nu_pi
                            - half_nu_plus_half * std::log(1.0 + z_sq_over_nu)
                            + log_const;
      acc += (long double)(z[i] * log_dens);
    }
    neg_ll -= (double)acc;
  }
  return neg_ll;
}

// ===========================================================================
// ANALYTIC GRADIENT OF THE CM-2 OBJECTIVE      (added 2026-08-31)
//
// CM-2 minimises, for component k,
//     F(gamma) = - sum_j sum_i z_i log f_FS(d_ij ; nu_j, gamma_j)
// with the Fernandez-Steel density written exactly as obj_gamma_cpp writes it:
//     zz   = d_ij * gamma_j^{-sign(d_ij)}          (sign(0) taken as +1)
//     logf = lgamma((nu+1)/2) - lgamma(nu/2) - log(nu pi)/2
//            - ((nu+1)/2) log(1 + zz^2/nu) + log 2 - log(gamma + 1/gamma).
//
// Only two terms depend on gamma. Using d zz / d gamma = -sign(d) zz / gamma,
//
//     dF/dgamma_j = -(1/gamma_j) sum_i z_i sign(d_ij) (nu_j+1) zz^2/(nu_j+zz^2)
//                   + n_z (1 - gamma_j^{-2}) / (gamma_j + gamma_j^{-1}),
//
// where n_z = sum_i z_i. Both sums run over the same (i, j) loop as the
// objective, so value and all m partial derivatives come out of ONE pass over
// the data. nlminb was differencing this objective, which costs m+1 passes per
// gradient; at m = 15 that is a sixteenfold difference on the step that
// profiling found to be 40-55% of total runtime.
//
// The objective is C^1 in gamma: at d_ij = 0 we have zz = 0 for every gamma,
// the log term contributes log(1) = 0, and its derivative vanishes from both
// sides, so the value assigned to sign(0) does not matter.
// ===========================================================================

// [[Rcpp::export]]
List obj_grad_gamma_cpp(Rcpp::NumericMatrix D_, Rcpp::NumericVector z_,
                        Rcpp::NumericVector nu_, Rcpp::NumericVector gamma_) {
  const int m = D_.nrow();
  const int N = D_.ncol();
  const double* D = D_.begin();          // column-major, D(j,i) = D[j + i*m]
  const double* z = z_.begin();

  long double nz = 0.0L;
  for (int i = 0; i < N; ++i) nz += (long double)z[i];
  const double n_z = (double)nz;

  double neg_ll = 0.0;
  Rcpp::NumericVector grad(m);

  for (int j = 0; j < m; ++j) {
    const double nu = nu_[j];
    const double gamma = gamma_[j];
    const double inv_gamma = 1.0 / gamma;
    const double log_const = LOG_2 - std::log(gamma + inv_gamma);
    const double half_nu = nu * 0.5;
    const double half_nu_plus_half = (nu + 1.0) * 0.5;
    const double lgamma_diff = fast_lgamma(half_nu_plus_half) - fast_lgamma(half_nu);
    const double log_nu_pi = 0.5 * std::log(nu * M_PI);

    long double acc = 0.0L;      // the objective
    long double gacc = 0.0L;     // sum_i z_i sign(d) (nu+1) zz^2 / (nu + zz^2)
    for (int i = 0; i < N; ++i) {
      const double d_val = D[j + (std::size_t)i * m];
      const bool pos = (d_val >= 0);
      const double scaling = pos ? inv_gamma : gamma;
      const double zz = d_val * scaling;
      const double zz2 = zz * zz;
      const double log_dens = lgamma_diff - log_nu_pi
                            - half_nu_plus_half * std::log(1.0 + zz2 / nu)
                            + log_const;
      acc += (long double)(z[i] * log_dens);
      const double sgn = pos ? 1.0 : -1.0;
      gacc += (long double)(z[i] * sgn * (nu + 1.0) * zz2 / (nu + zz2));
    }
    neg_ll -= (double)acc;
    grad[j] = -(double)gacc / gamma
              + n_z * (1.0 - inv_gamma * inv_gamma) / (gamma + inv_gamma);
  }

  return Rcpp::List::create(Rcpp::Named("value") = neg_ll,
                            Rcpp::Named("grad")  = grad);
}

// ===========================================================================
// COORDINATEWISE CM-2                                     (added 2026-09-01)
//
// The CM-2 objective is a sum over latent axes,
//     F(gamma) = sum_j F_j(gamma_j),
//     F_j(g) = - sum_i z_i [ log 2 - log(g + 1/g) + log t_nu_j( d_ij g^{-sign(d_ij)} ) ],
// so gamma_j appears in F_j alone and the m-dimensional block is m separate
// scalar problems. Solving it jointly is solving a harder problem than the
// model poses: measured on matched sub-problems, a coordinatewise search
// reached a lower objective than the joint quasi-Newton in 144 of 144 cases,
// and end to end it estimates gamma significantly more accurately
// (paired RMSE -0.0146, p = 0.001 over 120 fits). It does not change the
// partitions.
//
// Each coordinate is searched on its own profile in theta = log(gamma - 1):
//   1. a coarse grid over [lo, hi], which is a GLOBAL search of that profile
//      and does not depend on where the coordinate started;
//   2. bisection on F_j'(theta) inside the winning grid cell, when the
//      derivative brackets a stationary point there;
//   3. the incoming value is kept as a candidate, so each coordinate is
//      non-worsening by construction -- separability makes the guard exact
//      per coordinate rather than per vector.
//
// The grid must reach far to the left: gamma -> 1 is an exactly symmetric
// axis and sends theta -> -infinity. At convergence roughly half the fitted
// theta lie below -6 and the smallest observed was -16, so a range that stops
// at -6 silently clamps them.
// ===========================================================================

// value and d/dtheta of F_j at theta, for one latent axis
static inline void cm2_axis(const double* D, int m, int N, int j,
                            const double* z, double n_z, double nu, double theta,
                            double* val, double* der) {
  const double gamma = 1.0 + std::exp(theta);
  const double inv_gamma = 1.0 / gamma;
  const double log_const = LOG_2 - std::log(gamma + inv_gamma);
  const double half_nu_plus_half = (nu + 1.0) * 0.5;
  const double lgamma_diff = fast_lgamma(half_nu_plus_half) - fast_lgamma(nu * 0.5);
  const double log_nu_pi = 0.5 * std::log(nu * M_PI);

  long double acc = 0.0L, gacc = 0.0L;
  for (int i = 0; i < N; ++i) {
    const double d_val = D[j + (std::size_t)i * m];
    const bool pos = (d_val >= 0);
    const double zz = d_val * (pos ? inv_gamma : gamma);
    const double zz2 = zz * zz;
    acc += (long double)(z[i] * (lgamma_diff - log_nu_pi
                                 - half_nu_plus_half * std::log(1.0 + zz2 / nu)
                                 + log_const));
    if (der != nullptr)
      gacc += (long double)(z[i] * (pos ? 1.0 : -1.0)
                            * (nu + 1.0) * zz2 / (nu + zz2));
  }
  *val = -(double)acc;
  if (der != nullptr) {
    const double dF_dgamma = -(double)gacc / gamma
      + n_z * (1.0 - inv_gamma * inv_gamma) / (gamma + inv_gamma);
    *der = dF_dgamma * (gamma - 1.0);          // d gamma / d theta = gamma - 1
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector cm2_coord_cpp(Rcpp::NumericMatrix D_, Rcpp::NumericVector z_,
                                  Rcpp::NumericVector nu_,
                                  Rcpp::NumericVector theta0_,
                                  double lo = -20.0, double hi = 4.0,
                                  int ngrid = 49, int nbisect = 60) {
  const int m = D_.nrow(), N = D_.ncol();
  const double* D = D_.begin();
  const double* z = z_.begin();

  long double nzl = 0.0L;
  for (int i = 0; i < N; ++i) nzl += (long double)z[i];
  const double n_z = (double)nzl;

  Rcpp::NumericVector out(m);
  double v, d;

  for (int j = 0; j < m; ++j) {
    const double nu = nu_[j];

    // the incoming value is always a candidate
    double best_t = theta0_[j];
    cm2_axis(D, m, N, j, z, n_z, nu, best_t, &v, nullptr);
    double best_v = v;

    // 1. global scan of this axis's own profile
    int bi = -1; double bg = R_PosInf;
    const double step = (hi - lo) / (ngrid - 1);
    for (int g = 0; g < ngrid; ++g) {
      const double t = lo + g * step;
      cm2_axis(D, m, N, j, z, n_z, nu, t, &v, nullptr);
      if (v < bg) { bg = v; bi = g; }
    }
    if (bg < best_v) { best_v = bg; best_t = lo + bi * step; }

    // 2. bisection on the derivative inside the winning cell
    double a = lo + std::max(0, bi - 1) * step;
    double b = lo + std::min(ngrid - 1, bi + 1) * step;
    double da, db;
    cm2_axis(D, m, N, j, z, n_z, nu, a, &v, &da);
    cm2_axis(D, m, N, j, z, n_z, nu, b, &v, &db);
    if (da < 0.0 && db > 0.0) {                 // a stationary point is bracketed
      for (int it = 0; it < nbisect; ++it) {
        const double mid = 0.5 * (a + b);
        cm2_axis(D, m, N, j, z, n_z, nu, mid, &v, &d);
        if (d > 0.0) b = mid; else a = mid;
        if (b - a < 1e-12) break;
      }
      const double t = 0.5 * (a + b);
      cm2_axis(D, m, N, j, z, n_z, nu, t, &v, nullptr);
      if (v < best_v) { best_v = v; best_t = t; }
    }

    out[j] = best_t;
  }
  return out;
}

// ===========================================================================
// ANALYTIC GRADIENT OF THE CM-4 OBJECTIVE      (added 2026-08-30)
//
// CM-4 minimises, for component k,
//     J = n_k log|U| + 1/2 sum_ij a_ij Dt_ij^2 + lam_O sum_{r<s} S_rs^2
//                                              + lam_U sum_{r<s} U_rs^2
// with  V = U^{-1},  M = V O',  Yc = Y - 1 mu',  Dt = Yc M   (N x m), and
//     a_ij = z_i w_ij gamma_j^{-2 sign(Dt_ij)}.
//
// Writing G = A o Dt (elementwise) and G_M = Yc' G,
//     dJ/dmu      = -M colSums(G)
//     dJ/dU |wss  = -U^{-T} G_M M^T        (+ n_k on the log-diagonal params,
//                                            since n_k log|U| = n_k sum_j u_j)
//     G_O         = G_M' U^{-1}
//     dJ/dS       = -(I+P)' (R G_O) T',    T = (I+S)^{-1}, P = (I-S)T, O = R P
//     dJ/ds_rs    = (dJ/dS)_rs - (dJ/dS)_sr
// For the shared-structure modes the shared blocks accumulate over components:
//   EE ("full"): G_M = sum_k Yc_k' G_k, log-determinant contributes sum_k n_k
//   EV ("U")   : dJ/dU = sum_k -U^{-T} G_Mk M_k^T ; each O_k gets its own chain
//
// sign(Dt_ij) is evaluated at the CURRENT Dt on every call, never frozen at the
// incoming value: the exponent is part of e_ij^2 and moves with the parameters,
// so freezing it would optimise a surrogate instead of Q and forfeit the ascent
// guarantee. J is C^1 across {Dt_ij = 0}: a cell contributes c Dt^2 with
// c = gamma^-2 on one side and gamma^2 on the other, and both one-sided
// derivatives vanish at 0, so the gradient below is exact everywhere and the
// value assigned to sign(0) is immaterial. What jumps there is the SECOND
// derivative, so a quasi-Newton step taken near such a point can still be poor.
// The R driver therefore
// wraps every gradient solve in a monotonicity guard: if the returned point does
// not improve J, it falls back to the derivative-free optimizer, and failing
// that keeps the incoming parameters. ECM's monotone ascent of ell_p is then
// guaranteed regardless of what the optimizer did. Verified: the guard did not
// fire once in 24 complete fits.
//
// Verification:
//   gradient vs numDeriv           1e-9 .. 1e-5 relative, all three modes,
//                                  m in {2,3,5,8}, both reflection branches;
//                                  the residual is numDeriv's own O(d)
//                                  truncation error (error scales with the step)
//   value vs the existing objective  <= 1e-15 relative
//   CM-4 subproblem from a common start: same optimum, 1e-9 .. 7e-6 relative
//   end-to-end: ell_p differences have no direction (13 better / 11 worse,
//               p = 0.888); dARI -0.006, p = 0.347
// ===========================================================================

static arma::mat skew_from_upper(const arma::vec& s, int m) {
  arma::mat S(m, m, arma::fill::zeros);
  int idx = 0;
  for (int j = 0; j < m; ++j)
    for (int i = 0; i < j; ++i)
      if (idx < (int)s.n_elem) S(i, j) = s(idx++);
  return S - S.t();
}

// dJ/ds for one Cayley block, given dJ/dO; writes into out[offset ...]
static void cayley_chain(const arma::mat& G_O, const arma::mat& S, bool refl,
                         double lambda_O, const arma::vec& s_pars,
                         arma::vec& out, int offset) {
  const int m = S.n_rows;
  arma::mat I = arma::eye(m, m);
  arma::mat T = arma::inv(I + S);
  arma::mat P = (I - S) * T;
  arma::mat X = G_O; if (refl) X.row(0) *= -1.0;
  arma::mat gS = -((I + P).t() * X * T.t());
  int idx = 0;
  for (int j = 0; j < m; ++j)
    for (int i = 0; i < j; ++i) {
      out(offset + idx) = (gS(i, j) - gS(j, i)) + 2.0 * lambda_O * s_pars(idx);
      ++idx;
    }
}

// fills G = A o Dt and returns the weighted sum of squares for one component
static double accumulate_G(const arma::mat& Dt, const arma::vec& z,
                           const arma::mat& w, const arma::vec& g_sq,
                           const arma::vec& inv_g_sq, arma::mat& G) {
  const int N = Dt.n_rows, m = Dt.n_cols;
  double wss = 0.0;
  for (int j = 0; j < m; ++j) {
    const double* dp = Dt.colptr(j);
    const double* wp = w.colptr(j);
    const double* zp = z.memptr();
    double* gp = G.colptr(j);
    double gj = g_sq(j), igj = inv_g_sq(j);
    for (int i = 0; i < N; ++i) {
      double a = zp[i] * wp[i] * ((dp[i] >= 0) ? igj : gj);
      gp[i] = a * dp[i];
      wss += a * dp[i] * dp[i];
    }
  }
  return wss;
}

static inline double reg_U_upper(const arma::mat& U, int m, double lambda_U) {
  double r = 0.0;
  for (int j = 1; j < m; ++j)
    for (int i = 0; i < j; ++i) r += lambda_U * U(i, j) * U(i, j);
  return r;
}

// [[Rcpp::export]]
List obj_grad_hetero_cpp(arma::vec par, const arma::mat& Y, arma::vec gamma,
                         const arma::vec& z_k, const arma::mat& w,
                         double lambda_U, double lambda_O, bool use_reflection) {
  const int m = Y.n_cols, N = Y.n_rows;
  const int nA = (m > 1) ? m * (m - 1) / 2 : 0;
  const double n_z = arma::sum(z_k);

  arma::vec mu = par.subvec(0, m - 1);
  arma::vec u_pars = par.subvec(m, 2 * m + nA - 1);
  arma::vec s_pars = (nA > 0) ? par.subvec(2 * m + nA, 2 * m + 2 * nA - 1)
                              : arma::vec();

  arma::mat U = build_U_cpp(u_pars, m);
  if (!is_well_conditioned(U))
    return List::create(Named("value") = PENALTY_VALUE,
                        Named("grad") = arma::zeros(par.n_elem));

  arma::mat S = (m > 1) ? skew_from_upper(s_pars, m) : arma::zeros(m, m);
  arma::mat I = arma::eye(m, m);
  arma::mat P = (I - S) * arma::inv(I + S);
  arma::mat O = P; if (use_reflection) O.row(0) *= -1.0;

  arma::mat Vi;
  if (!arma::inv(Vi, arma::trimatu(U)))
    return List::create(Named("value") = PENALTY_VALUE,
                        Named("grad") = arma::zeros(par.n_elem));
  arma::mat M = Vi * O.t();
  arma::mat Yc = Y.each_row() - mu.t();
  arma::mat Dt = Yc * M;

  arma::vec g_sq = gamma % gamma, inv_g_sq = 1.0 / g_sq;
  arma::mat G(N, m);
  double wss = accumulate_G(Dt, z_k, w, g_sq, inv_g_sq, G);

  double log_det_U = arma::sum(arma::log(U.diag()));
  double reg = 0.0;
  if (m > 1) { reg = reg_U_upper(U, m, lambda_U);
               reg += lambda_O * arma::dot(s_pars, s_pars); }
  double value = n_z * log_det_U + 0.5 * wss + reg;

  arma::vec grad(par.n_elem, arma::fill::zeros);
  arma::rowvec cs = arma::sum(G, 0);
  grad.subvec(0, m - 1) = -(M * cs.t());
  arma::mat G_M = Yc.t() * G;
  arma::mat gU = -(Vi.t() * G_M * M.t());
  for (int i = 0; i < m; ++i) grad(m + i) = gU(i, i) * U(i, i) + n_z;
  { int idx = m;
    for (int j = 0; j < m; ++j)
      for (int i = 0; i < j; ++i) {
        grad(m + idx) = gU(i, j) + 2.0 * lambda_U * U(i, j); ++idx; } }
  if (m > 1) cayley_chain(G_M.t() * Vi, S, use_reflection, lambda_O, s_pars,
                          grad, 2 * m + nA);
  return List::create(Named("value") = value, Named("grad") = grad);
}

// [[Rcpp::export]]
List obj_grad_full_cpp(arma::vec par, const arma::mat& Y, List components,
                       const arma::mat& Z, const arma::cube& W,
                       double lambda_U, double lambda_O, bool use_reflection) {
  const int m = Y.n_cols, N = Y.n_rows, K = components.size();
  const int nA = (m > 1) ? m * (m - 1) / 2 : 0;
  int off = K * m;

  arma::vec u_pars = par.subvec(off, off + m + nA - 1);
  arma::vec s_pars = (nA > 0) ? par.subvec(off + m + nA, off + m + 2 * nA - 1)
                              : arma::vec();
  arma::mat U = build_U_cpp(u_pars, m);
  if (!is_well_conditioned(U))
    return List::create(Named("value") = PENALTY_VALUE,
                        Named("grad") = arma::zeros(par.n_elem));
  arma::mat S = (m > 1) ? skew_from_upper(s_pars, m) : arma::zeros(m, m);
  arma::mat I = arma::eye(m, m);
  arma::mat P = (I - S) * arma::inv(I + S);
  arma::mat O = P; if (use_reflection) O.row(0) *= -1.0;

  arma::mat Vi;
  if (!arma::inv(Vi, arma::trimatu(U)))
    return List::create(Named("value") = PENALTY_VALUE,
                        Named("grad") = arma::zeros(par.n_elem));
  arma::mat M = Vi * O.t();
  double log_det_U = arma::sum(arma::log(U.diag()));

  arma::vec grad(par.n_elem, arma::fill::zeros);
  arma::mat G_M(m, m, arma::fill::zeros);
  double total = 0.0, n_tot = 0.0;

  for (int k = 0; k < K; ++k) {
    arma::vec z = Z.col(k);
    double n_k = arma::sum(z);
    if (n_k < 1e-8) continue;
    n_tot += n_k;
    List th = components[k];
    arma::vec gam = as<arma::vec>(th["gamma"]);
    arma::vec g_sq = gam % gam, inv_g_sq = 1.0 / g_sq;
    arma::vec mu = par.subvec(k * m, k * m + m - 1);
    arma::mat Yc = Y.each_row() - mu.t();
    arma::mat Dt = Yc * M;
    arma::mat G(N, m);
    double wss = accumulate_G(Dt, z, W.slice(k), g_sq, inv_g_sq, G);
    total += n_k * log_det_U + 0.5 * wss;
    arma::rowvec cs = arma::sum(G, 0);
    grad.subvec(k * m, k * m + m - 1) = -(M * cs.t());
    G_M += Yc.t() * G;
  }

  double reg = 0.0;
  if (m > 1) { reg = reg_U_upper(U, m, lambda_U);
               reg += lambda_O * arma::dot(s_pars, s_pars); }

  arma::mat gU = -(Vi.t() * G_M * M.t());
  for (int i = 0; i < m; ++i) grad(off + i) = gU(i, i) * U(i, i) + n_tot;
  { int idx = m;
    for (int j = 0; j < m; ++j)
      for (int i = 0; i < j; ++i) {
        grad(off + idx) = gU(i, j) + 2.0 * lambda_U * U(i, j); ++idx; } }
  if (m > 1) cayley_chain(G_M.t() * Vi, S, use_reflection, lambda_O, s_pars,
                          grad, off + m + nA);
  return List::create(Named("value") = total + reg, Named("grad") = grad);
}

// [[Rcpp::export]]
List obj_grad_u_only_cpp(arma::vec par, const arma::mat& Y, List components,
                         const arma::mat& Z, const arma::cube& W,
                         double lambda_U, double lambda_O,
                         LogicalVector use_reflection_vec) {
  const int m = Y.n_cols, N = Y.n_rows, K = components.size();
  const int nA = (m > 1) ? m * (m - 1) / 2 : 0;
  int uoff = K * m + K * nA;

  arma::vec u_pars = par.subvec(uoff, uoff + m + nA - 1);
  arma::mat U = build_U_cpp(u_pars, m);
  if (!is_well_conditioned(U))
    return List::create(Named("value") = PENALTY_VALUE,
                        Named("grad") = arma::zeros(par.n_elem));
  arma::mat Vi;
  if (!arma::inv(Vi, arma::trimatu(U)))
    return List::create(Named("value") = PENALTY_VALUE,
                        Named("grad") = arma::zeros(par.n_elem));
  double log_det_U = arma::sum(arma::log(U.diag()));

  arma::vec grad(par.n_elem, arma::fill::zeros);
  arma::mat gU_acc(m, m, arma::fill::zeros);
  arma::mat I = arma::eye(m, m);
  double total = 0.0, n_tot = 0.0, regO = 0.0;

  for (int k = 0; k < K; ++k) {
    arma::vec z = Z.col(k);
    double n_k = arma::sum(z);
    if (n_k < 1e-8) continue;
    n_tot += n_k;
    arma::vec s_k = (nA > 0) ? par.subvec(K * m + k * nA, K * m + k * nA + nA - 1)
                             : arma::vec();
    arma::mat S = (m > 1) ? skew_from_upper(s_k, m) : arma::zeros(m, m);
    arma::mat P = (I - S) * arma::inv(I + S);
    arma::mat O = P; if (use_reflection_vec[k]) O.row(0) *= -1.0;
    arma::mat M = Vi * O.t();

    List th = components[k];
    arma::vec gam = as<arma::vec>(th["gamma"]);
    arma::vec g_sq = gam % gam, inv_g_sq = 1.0 / g_sq;
    arma::vec mu = par.subvec(k * m, k * m + m - 1);
    arma::mat Yc = Y.each_row() - mu.t();
    arma::mat Dt = Yc * M;
    arma::mat G(N, m);
    double wss = accumulate_G(Dt, z, W.slice(k), g_sq, inv_g_sq, G);
    if (m > 1) regO += lambda_O * arma::dot(s_k, s_k);
    total += n_k * log_det_U + 0.5 * wss;

    arma::rowvec cs = arma::sum(G, 0);
    grad.subvec(k * m, k * m + m - 1) = -(M * cs.t());
    arma::mat G_Mk = Yc.t() * G;
    gU_acc += -(Vi.t() * G_Mk * M.t());
    if (m > 1) cayley_chain(G_Mk.t() * Vi, S, use_reflection_vec[k], lambda_O,
                            s_k, grad, K * m + k * nA);
  }

  double regU = (m > 1) ? reg_U_upper(U, m, lambda_U) : 0.0;
  for (int i = 0; i < m; ++i) grad(uoff + i) = gU_acc(i, i) * U(i, i) + n_tot;
  { int idx = m;
    for (int j = 0; j < m; ++j)
      for (int i = 0; i < j; ++i) {
        grad(uoff + idx) = gU_acc(i, j) + 2.0 * lambda_U * U(i, j); ++idx; } }
  return List::create(Named("value") = total + regO + regU, Named("grad") = grad);
}
