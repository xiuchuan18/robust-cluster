# Section6_ScaleLib.R
#
# Data and timing helpers for the Section 6.3 scalability study.

gen_scale <- function(N, m, K, seed, sep = 5) {
  set.seed(seed)
  h <- max(1, min(3, floor(m / 2)))
  mu <- lapply(seq_len(K), function(k)
    if (k == 1) rep(0, m) else mu_at(m, ((k - 2) * h + 1):((k - 1) * h), sep))
  gm <- lapply(seq_len(K), function(k)
    seq(0.55 + 0.2 * (k - 1), 2.0 - 0.15 * (k - 1), length.out = m))
  nu <- lapply(seq_len(K), function(k) seq(3 + k, 6 + 2 * k, length.out = m))
  U  <- diag(exp(seq(-0.25, 0.25, length.out = m)))
  if (m > 1) U[upper.tri(U)] <- runif(m * (m - 1) / 2, -0.25, 0.25)
  O  <- lapply(seq_len(K), function(k) qr.Q(qr(matrix(rnorm(m^2), m, m))))
  nk <- as.vector(rmultinom(1, N, rep(1 / K, K)))
  Y <- do.call(rbind, lapply(seq_len(K), function(k)
    assemble(nk[k], mu[[k]], U, O[[k]], nu[[k]], gm[[k]])))
  p <- sample.int(N)
  list(Y = Y[p, , drop = FALSE], labels = rep(seq_len(K), nk)[p])
}

# One timed fit. For the MixSIT variants the base fit and the flip polish are
# timed separately and the ECM iteration count is taken from the fitter itself,
# so cost per iteration is measured rather than inferred.
# Competitors go through fit_multistart() so Section 6.3 uses the same protocol
# as Sections 4 and 5.
#
# MixSIT keeps its own path because 6.3 needs the ECM time and the flip time
# separately, which is a decomposition fit_multistart does not expose.
time_fit <- function(method, Y, K, lam, tag = NULL, rep = NULL) {
  d <- ncol(Y); n <- nrow(Y)
  if (!grepl("MixSIT", method)) {
    t0 <- Sys.time()
    r <- tryCatch(fit_multistart(method, Y, K, R_RESTART, COMPETITOR_MAX_ITER,
                                 tag = tag, rep = rep), error = function(e) NULL)
    tb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    it <- if (is.null(r) || is.null(r$fit)) NA_integer_ else
          tryCatch(switch(method,
            MSGHD = length(r$fit@loglik),
            GMM   = NA_integer_,
            NA_integer_), error = function(e) NA_integer_)
    return(list(ok = !is.null(r) && !is.null(r$cluster), t_base = tb,
                t_flip = 0, iters = it, flips = NA_integer_))
  }
  mode <- c(MixSIT_None = "none", MixSIT_U = "U", MixSIT_Full = "full")[[method]]
  t0 <- Sys.time()
  f0 <- tryCatch(fitMixSIT_ECM(Y, K = K, constraint_mode = mode,
          n_short = 1, short_iter = 1, init_methods = "kmeans",
          max_iter = 500, tol = MIXSIT$tol, lambda_U = 0, lambda_O = 0,
          lambda_magnitude = 0, lambda_hetero = lam, flips = FALSE), error = function(e) NULL)
  tb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(f0)) return(list(ok = FALSE, t_base = tb, t_flip = 0,
                               iters = NA_integer_, flips = NA_integer_))
  t0 <- Sys.time()
  f1 <- tryCatch(polish_from(Y, f0$params, K, lam, mode), error = function(e) NULL)
  tf <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(ok = !is.null(f1), t_base = tb, t_flip = tf,
       iters = f0$iterations, flips = if (is.null(f1)) NA_integer_ else f1$n_flips)
}
