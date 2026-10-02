# Section4_Lib.R
#
# Final generators, information criteria, and fitting protocol for the six Section 4 scenarios. Each replicate is run in an isolated R process.

suppressPackageStartupMessages({
  library(mvtnorm); library(mclust); library(EMMIXskew); library(EMMIXuskew)
  library(MixGHD); library(ContaminatedMixt); library(MASS)
})

# Every method is fitted once. Non-GMM methods use one package-native k-means initializer with nstart=5; GMM uses mclust's deterministic default.
R_RESTART <- 1

CONFIG_SIM <- list(
  n_sim = 100, K = 3, seed = SEED_SIM, max_iter = MIXSIT$max_iter,
  scenarios = list(
    # Internal names are stable identifiers; displayed labels are below.
    S1_Gaussian     = list(N = 600, d = 5),
    S2_Imbalanced   = list(N = 600, d = 5),
    # S3: heterogeneous tail indices on the OBSERVED coordinates, weakly dependent (Gaussian copula, equicorrelation 0.30, t margins with S3's
    # per-axis nu). Replaces the t-copula cell, whose heaviness came from one scale shared across a row -- tMM's own model -- and which MixSIT lost.
    S3_GaussCopula  = list(N = 600, d = 5),
    S4_Skewed       = list(N = 600, d = 5),
    # S5 (manuscript): cell-level contamination at m = 10, unchanged generator.
    S7_CellWise10   = list(N = 700, d = 10),
    # S6: aligned geometry with wider per-axis skewness (gamma^1.5) and
    # moderate tails (nu in [5, 8] on every axis).
    S6_AlignedWide  = list(N = 600, d = 10)
  )
)

CELL_LABEL <- c(S1_Gaussian = "S1", S2_Imbalanced = "S2", S3_GaussCopula = "S3: Heavy tailed",
                S4_Skewed = "S4", S7_CellWise10 = "S5", S6_AlignedWide = "S6")

# Data generators

# Fernandez-Steel skew-t source, and mean-matched assembly
rfs <- function(n, nu, gamma) {
  pos <- runif(n) < gamma^2 / (1 + gamma^2)
  z <- numeric(n)
  z[pos]  <-  abs(rt(sum(pos),  df = nu)) * gamma
  z[!pos] <- -abs(rt(sum(!pos), df = nu)) / gamma
  z
}
fs_mean <- function(g, nu) {
  if (nu <= 1) return(0)
  2 * sqrt(nu) * gamma((nu + 1) / 2) / (sqrt(pi) * (nu - 1) * gamma(nu / 2)) * (g - 1/g)
}
# y = mode_k + U'O' eps, with mode_k chosen so the component MEAN is mu_target
assemble <- function(nk, mu_target, U, O, nu, gm) {
  m <- length(mu_target)
  E <- vapply(seq_len(m), function(j) rfs(nk, nu[j], gm[j]), numeric(nk))
  if (is.null(dim(E))) E <- matrix(E, nrow = nk)
  me <- vapply(seq_len(m), function(j) fs_mean(gm[j], nu[j]), numeric(1))
  mode_k <- mu_target - as.vector(t(U) %*% t(O) %*% me)
  sweep(t(t(U) %*% t(O) %*% t(E)), 2, mode_k, "+")
}
mu_at  <- function(d, pos, v) { x <- numeric(d); p <- pos[pos <= d]; x[p] <- v; x }

finish <- function(Y, lab, mask = NULL, seed_perm = NULL) {
  if (!is.null(seed_perm)) set.seed(seed_perm)
  p <- sample.int(nrow(Y))
  list(Y = Y[p, , drop = FALSE], labels = lab[p],
       contam = if (is.null(mask)) NULL else mask[p, , drop = FALSE]) }

generate_S1_gaussian <- function(N, d, K = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_k <- as.vector(rmultinom(1, N, rep(1/K, K)))
  mu_list <- list(rep(0, d), c(4, 4, rep(0, d - 2)), c(4, -2, rep(0, d - 2)))
  Y <- matrix(0, N, d); labels <- integer(N); idx <- 1
  for (k in 1:K) if (n_k[k] > 0) {
    Y[idx:(idx + n_k[k] - 1), ] <- mvtnorm::rmvnorm(n_k[k], mu_list[[k]], diag(d))
    labels[idx:(idx + n_k[k] - 1)] <- k; idx <- idx + n_k[k] }
  finish(Y, labels)
}

generate_S2_imbalanced <- function(N, d, K = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_k <- as.vector(rmultinom(1, N, c(0.1, 0.3, 0.6)))
  mu_list <- list(rep(0, d), c(2, 2, rep(0, d - 2)), c(-2, -2, rep(0, d - 2)))
  Y <- matrix(0, N, d); labels <- integer(N); idx <- 1
  for (k in 1:K) if (n_k[k] > 0) {
    Y[idx:(idx + n_k[k] - 1), ] <- mvtnorm::rmvnorm(n_k[k], mu_list[[k]], diag(d))
    labels[idx:(idx + n_k[k] - 1)] <- k; idx <- idx + n_k[k] }
  finish(Y, labels)
}

get_fixed_covariance <- function(d, type = 1) {
  if (type == 1) {
    rho <- 0.7; Sigma <- rho^abs(outer(1:d, 1:d, "-"))
    dv <- seq(0.5, 1.5, length.out = d)
    diag(sqrt(dv)) %*% Sigma %*% diag(sqrt(dv))
  } else if (type == 2) {
    Sigma <- matrix(0.5, d, d); diag(Sigma) <- 1
    dv <- seq(1, 2, length.out = d)
    diag(sqrt(dv)) %*% Sigma %*% diag(sqrt(dv))
  } else {
    rho <- -0.6; Sigma <- rho^abs(outer(1:d, 1:d, "-")); diag(Sigma) <- 1.5
    Sigma
  }
}

# correct sampler for the unrestricted skew-t of Sahu, Dey & Branco
rmst_sdb <- function(n, mu, sigma, delta, dof) {
  p <- length(mu)
  U1 <- abs(MASS::mvrnorm(n, rep(0, p), diag(p)))
  U0 <- MASS::mvrnorm(n, rep(0, p), sigma)
  if (p == 1) { U1 <- matrix(U1, ncol = 1); U0 <- matrix(U0, ncol = 1) }
  G <- stats::rgamma(n, dof / 2, dof / 2)
  sweep((sweep(U1, 2, delta, "*") + U0) / sqrt(G), 2, mu, "+")
}

S3_NU_AXIS <- list(c(3, 5, 8, 15, 30), c(30, 15, 8, 5, 3), c(8, 30, 3, 15, 5))
S3_CENTRE  <- list(c(0, 0), c(1, 1), c(-1, -1.5))
S3_SCALE   <- 1.50      # calibrated so the class separation matches the 2.30 of the published design; it was not chosen on any outcome

.gen_uskew <- function(N, d, K, seed, nu_val, mu_list, delta_list) {
  if (!is.null(seed)) set.seed(seed)
  n_k <- as.vector(rmultinom(1, N, rep(1/K, K)))
  sigma_list <- list(get_fixed_covariance(d, 1), get_fixed_covariance(d, 2),
                     get_fixed_covariance(d, 3))
  Y <- matrix(0, N, d); labels <- integer(N); idx <- 1
  for (k in 1:K) if (n_k[k] > 0) {
    Y[idx:(idx + n_k[k] - 1), ] <- rmst_sdb(
      n = n_k[k], mu = mu_list[[k]], sigma = sigma_list[[k]],
      delta = delta_list[[k]], dof = nu_val)
    labels[idx:(idx + n_k[k] - 1)] <- k; idx <- idx + n_k[k] }
  finish(Y, labels)
}
# S3: Gaussian copula, t margins with heterogeneous nu
# z ~ N(0, R) with R equicorrelated at S3_GC_RHO; u_j = Phi(z_j); x_j = F^-1_{t_nu_j}(u_j) / sqrt(nu_j / (nu_j - 2)), unit-variance margins.
# The heaviness lives on the observed coordinates -- no latent-axis construction -- and each coordinate has its own tail index, S3_NU_AXIS as
# before. Centres and scale as in the t-copula cell. Density is closed form (Section4_BayesCeiling.R, logdens_gc).
S3_GC_RHO <- 0.30
generate_S3_gausscopula <- function(N, d = 5, K = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_k <- as.vector(rmultinom(1, N, rep(1 / K, K)))
  R <- matrix(S3_GC_RHO, d, d); diag(R) <- 1; L <- chol(R)
  Y <- matrix(0, N, d); labels <- integer(N); idx <- 1
  for (k in seq_len(K)) if (n_k[k] > 0) {
    nk <- S3_NU_AXIS[[k]][1:d]; n <- n_k[k]
    U <- stats::pnorm(matrix(stats::rnorm(n * d), n, d) %*% L)
    U <- pmin(pmax(U, 1e-12), 1 - 1e-12)
    X <- vapply(seq_len(d), function(j) stats::qt(U[, j], df = nk[j]) / sqrt(nk[j] / (nk[j] - 2)), numeric(n))
    mu <- S3_SCALE * mu_at(d, 1:2, S3_CENTRE[[k]])
    Y[idx:(idx + n - 1), ] <- sweep(X, 2, mu, "+")
    labels[idx:(idx + n - 1)] <- k; idx <- idx + n }
  finish(Y, labels)
}
generate_S4_skewed <- function(N, d = 5, K = 3, seed = NULL)
  .gen_uskew(N, d, K, seed, 20,
             list(rep(0, d), c(2, 2, rep(0, d - 2)), c(-2, -2, rep(0, d - 2))),
             list(c(10,6,5,-8,4)[1:d], c(12,4,8,-6,6)[1:d], c(10,-5,-4,-4,10)[1:d]))

# S6: aligned geometry, wider skewness, moderate tails
# Per-axis gamma is raised to the 1.5 power (log-spread x 1.5) and nu set to seq(5, 8) on every axis.
# Direction-specific skewness is varied while nu is held away from the shared-scale heavy-tail regime.
S6_NU <- function(d) rep(list(seq(5, 8, length.out = d)), 3)
S6_GM <- function(d) list(seq(0.6, 1.8, length.out = d)^1.5, seq(1.2, 2.2, length.out = d)^1.5,
                          seq(0.5, 1.5, length.out = d)^1.5)
generate_S6_alignedwide <- function(N, d = 10, K = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_k <- as.vector(rmultinom(1, N, rep(1/K, K)))
  mu_list <- list(rep(0, d), mu_at(d, 1:3, 2), mu_at(d, 4:6, 2))
  nu_list <- S6_NU(d); gm_list <- S6_GM(d)
  O_list <- lapply(1:K, function(k) qr.Q(qr(matrix(rnorm(d^2), d, d))))
  Y <- matrix(0, N, d); labels <- integer(N); idx <- 1
  for (k in 1:K) if (n_k[k] > 0) {
    Y[idx:(idx + n_k[k] - 1), ] <- assemble(n_k[k], mu_list[[k]], diag(d),
                                            O_list[[k]], nu_list[[k]], gm_list[[k]])
    labels[idx:(idx + n_k[k] - 1)] <- k; idx <- idx + n_k[k] }
  finish(Y, labels)
}

# S5: cell-wise contamination, m = 10
generate_cellwise_dim <- function(N, d, K = 3, alpha = 0.05, eta = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_k <- as.vector(rmultinom(1, N, rep(1/K, K)))
  h <- max(2, min(4, floor(d / 2)))
  mu_list <- list(numeric(d), mu_at(d, 1:h, 2.5), mu_at(d, (h + 1):(2 * h), 2.5))
  Y <- matrix(0, N, d); labels <- integer(N); idx <- 1
  for (k in 1:K) if (n_k[k] > 0) {
    Y[idx:(idx + n_k[k] - 1), ] <- mvtnorm::rmvnorm(n_k[k], mu_list[[k]], diag(d))
    labels[idx:(idx + n_k[k] - 1)] <- k; idx <- idx + n_k[k] }
  # contamination on its own stream, so alpha cannot move the clean sample
  if (!is.null(seed)) set.seed(mk_seed(seed, "contam"))
  C <- matrix(runif(N * d) < alpha, N, d)
  Y[C] <- Y[C] + rnorm(sum(C), 0, eta)
  finish(Y, labels, C, seed_perm = if (is.null(seed)) NULL else mk_seed(seed, "perm"))
}

generate_cell <- function(cell, N, d, K, seed) switch(cell,
  S1_Gaussian     = generate_S1_gaussian(N, d, K, seed),
  S2_Imbalanced   = generate_S2_imbalanced(N, d, K, seed),
  S3_GaussCopula  = generate_S3_gausscopula(N, d, K, seed),
  S6_AlignedWide  = generate_S6_alignedwide(N, d, K, seed),
  S4_Skewed       = generate_S4_skewed(N, d, K, seed),
  S7_CellWise10   = generate_cellwise_dim(N, d, K, seed = seed),
  stop("unknown cell: ", cell))

CONTAM_CELLS <- c("S7_CellWise10")

# Information criteria: actual MSGHD block count and shared soft-posterior ICL.
count_parameters <- function(model_name, K, d) {
  base_params <- (K - 1) + (K * d)
  cov_params  <- K * d * (d + 1) / 2
  if (model_name == "GMM")   return(base_params + cov_params)
  if (model_name == "tMM")   return(base_params + cov_params + K)
  if (model_name == "SNMM")  return(base_params + cov_params + K * d)
  if (model_name == "STMM")  return(base_params + cov_params + K * d + K)
  if (model_name == "MSGHD") return((K - 1) + K * (d * (d - 1) / 2 + 5 * d))
  if (model_name == "CGMM")  return(base_params + cov_params + K * 2)
  if (grepl("MixSIT", model_name)) {
    gamma_nu_params <- 2 * K * d
    if (model_name == "MixSIT_None") {
      o_params <- K * d * (d - 1) / 2; u_params <- K * d * (d + 1) / 2
    } else if (model_name == "MixSIT_U") {
      o_params <- K * d * (d - 1) / 2; u_params <- d * (d + 1) / 2
    } else {
      o_params <- d * (d - 1) / 2;     u_params <- d * (d + 1) / 2
    }
    return((K - 1) + (K * d) + u_params + o_params + gamma_nu_params)
  }
  NA
}

calc_loglik <- function(Y, fit, model_name) {
  tryCatch({
    if (model_name == "MSGHD") {
      ll <- fit@loglik
      if (is.null(ll) || length(ll) == 0) return(NA_real_)
      return(as.numeric(tail(ll, 1)))
    }
    # MixSIT: the RAW log-likelihood, not the penalized objective the algorithm
    # ascends. fit$loglik is the ell_p trajectory; fit$loglik_obs is ell.
    loglik <- if (grepl("MixSIT", model_name))
                (if (!is.null(fit$loglik_obs)) fit$loglik_obs else tail(fit$loglik, 1))
              else if (model_name == "CGMM")   fit$models[[1]]$loglik
              else                             fit$loglik
    if (is.null(loglik) || length(loglik) == 0) return(NA_real_)
    as.numeric(tail(loglik, 1))
  }, error = function(e) NA_real_)
}

# Shared soft-posterior entropy convention for every information criterion. The numerical floor matches the original Section 4 / Section 6 convention.
# K=1 has exactly zero entropy; missing or invalid posteriors stay missing.
posterior_entropy <- function(z, n = NULL) {
  if (is.null(z)) return(NA_real_)
  z <- as.matrix(z)
  if (!is.null(n) && nrow(z) != n) return(NA_real_)
  if (!nrow(z) || !ncol(z) || any(!is.finite(z)) || any(z < 0)) return(NA_real_)
  rs <- rowSums(z)
  if (any(rs <= 0)) return(NA_real_)
  if (ncol(z) == 1L) return(0)
  zz <- pmax(z / rs, 1e-10)
  as.numeric(-sum(zz * log(zz)))
}

model_posterior <- function(fit, model_name, n = NULL, K = NULL) {
  tryCatch({
    if (model_name == "MSGHD") {
      # MixGHD's S4 vector slot drops the column-major matrix dimensions.
      if (is.null(n)) n <- nrow(fit@data)
      z <- fit@z
      if (is.null(K)) K <- length(z) / n
      if (length(z) != n * K || K < 1 || K != as.integer(K)) return(NULL)
      return(matrix(as.numeric(z), nrow = n, ncol = K))
    }
    if (model_name == "GMM") fit$z
    else if (model_name %in% c("tMM", "SNMM", "STMM")) fit$tau
    else if (model_name == "CGMM") fit$models[[1]]$posterior
    else if (grepl("MixSIT", model_name)) fit$posterior
    else NULL
  }, error = function(e) NULL)
}

calc_information_criteria <- function(Y, fit, model_name, K, d) {
  n <- nrow(Y)
  loglik <- calc_loglik(Y, fit, model_name)
  if (is.null(loglik) || is.na(loglik) || !is.finite(loglik))
    return(list(AIC = NA_real_, BIC = NA_real_, ICL = NA_real_))
  n_params <- count_parameters(model_name, K, d)
  if (is.na(n_params)) return(list(AIC = NA_real_, BIC = NA_real_, ICL = NA_real_))

  aic <- -2 * loglik + 2 * n_params
  bic <- -2 * loglik + n_params * log(n)
  z <- model_posterior(fit, model_name, n, K)
  ent <- posterior_entropy(z, n)
  list(AIC = as.numeric(aic), BIC = as.numeric(bic),
       ICL = if (!is.na(ent)) as.numeric(bic + 2 * ent) else NA_real_)
}

# One fit = R independent restarts, selected on the method's OWN log-likelihood EmSkew, MSGHD and CNmixt print their whole fitted object to stdout. Over
# 1100 replicates x 9 methods that is gigabytes of console I/O, which on Windows costs more wall time than the fits. Swallow it.
.quiet <- function(expr) {
  f <- tempfile()
  sink(f, type = "output")
  on.exit({ sink(type = "output"); unlink(f) }, add = TRUE)
  suppressMessages(expr)
}

# EMMIXskew's native k-means branch already uses nstart = 5, but init.mix() silently switches to ten random partitions when that start is unusable. The
# final protocol forbids changing initialization families. Build the identical package-native initialization explicitly, once, and pass it to EmSkew().
.emmix_kmeans_init <- function(Y, K, distr, ncov = 3L, initloop = 20L,
                               nstart = KMEANS_NSTART) {
  cl <- if (K == 1L) rep.int(1L, nrow(Y)) else
    stats::kmeans(Y, centers = K, nstart = nstart)$cluster
  if (K > 1L && min(tabulate(cl, nbins = K)) < 10L) return(NULL)
  ans <- try(EMMIXskew:::initEmmix(as.matrix(Y), K, cl, distr, ncov,
                                    initloop), silent = TRUE)
  if (inherits(ans, "try-error") || length(ans) != 8L ||
      isTRUE(as.logical(ans$error))) return(NULL)
  ans
}

# MSGHD exposes method = "km" but not kmeans()'s nstart. This reproduces MixGHD 2.3.7's rmgparMS() path exactly (partition, rgparMS construction, then
# ten EMgrstepMS warm-up iterations), changing only nstart from 1 to 5.
.msghd_kmeans_gpar <- function(Y, K, warmup = 10L,
                               nstart = KMEANS_NSTART) {
  if (K == 1L) {
    centers <- matrix(colMeans(Y), nrow = 1L)
    cl <- rep.int(1L, nrow(Y))
  } else {
    km <- stats::kmeans(Y, centers = K, nstart = nstart)
    centers <- km$centers
    cl <- km$cluster
  }
  z <- MixGHD:::combinewk(weights = matrix(0, nrow(Y), K), label = cl)
  gp <- MixGHD:::rgparMS(g = K, data = Y, w = z, l = centers)
  for (j in seq_len(warmup))
    gp <- MixGHD:::EMgrstepMS(data = Y, gpar = gp, v = 1,
                              label = cl, w = z, it = j)
  gp
}

fit_multistart <- function(code, Y, K, R, max_iter, tag = NULL, rep = NULL) {
  N <- nrow(Y); best <- NULL; best_ll <- -Inf
  keep <- function(f, ll) if (is.finite(ll) && ll > best_ll) {
    best_ll <<- ll; best <<- f }
  if (!is.null(tag) && !is.null(rep)) seed_for(tag, "fit", code, rep)

  if (code == "GMM") {
    f0 <- tryCatch(mclust::Mclust(Y, G = K, modelNames = "VVV", verbose = FALSE),
                   error = function(e) NULL)
    if (!is.null(f0)) keep(f0, f0$loglik)
    for (r in seq_len(R - 1)) {
      f <- tryCatch(mclust::Mclust(Y, G = K, modelNames = "VVV", verbose = FALSE,
             initialization = list(subset = sample.int(N, round(N * 0.6)))),
             error = function(e) NULL)
      if (!is.null(f)) keep(f, f$loglik)
    }
    return(list(fit = best, cluster = if (is.null(best)) NULL else best$classification))
  }

  if (code %in% c("tMM", "SNMM", "STMM")) {
    # Candidate protocol: one package-native k-means initialization and no random or hierarchical alternatives.
    dist <- c(tMM = "mvt", SNMM = "msn", STMM = "mst")[[code]]
    it <- max(max_iter, COMPETITOR_DEFAULTS$EmSkew$itmax)
    ep <- min(COMPETITOR_TOL, COMPETITOR_DEFAULTS$EmSkew$epsilon)
    ini <- .quiet(tryCatch(.emmix_kmeans_init(
      Y, K, dist, ncov = 3L,
      initloop = COMPETITOR_DEFAULTS$EmSkew$initloop),
      error = function(e) NULL))
    if (is.null(ini)) return(list(fit = NULL, cluster = NULL))
    f <- .quiet(tryCatch(EMMIXskew::EmSkew(
        Y, g = K, distr = dist, ncov = 3, init = ini,
        itmax = it, epsilon = ep, debug = FALSE), error = function(e) NULL))
    return(list(fit = f, cluster = if (is.null(f)) NULL else f$clust))
  }

  if (code == "MSGHD") {
    gp <- .quiet(tryCatch(.msghd_kmeans_gpar(Y, K), error = function(e) NULL))
    if (is.null(gp)) return(list(fit = NULL, cluster = NULL))
    f <- .quiet(tryCatch(MixGHD::MSGHD(data = Y, gpar0 = gp, G = K,
           method = "km",
           max.iter = COMPETITOR_DEFAULTS$MSGHD$max.iter,
           eps = min(COMPETITOR_TOL, COMPETITOR_DEFAULTS$MSGHD$eps),
           scale = FALSE, modelSel = "ICL"),
           error = function(e) NULL))
    ll <- if (is.null(f)) -Inf else
          tryCatch(as.numeric(tail(f@loglik, 1)), error = function(e) -Inf)
    if (!is.null(f) && is.finite(ll)) keep(f, ll)
    return(list(fit = best, cluster = if (is.null(best)) NULL else best@map))
  }

  if (code == "CGMM") {
    for (ini in head(c("kmeans", "mixt", "random.post", "random.clas",
                       "random.post"), R)) {
      f <- .quiet(tryCatch(ContaminatedMixt::CNmixt(Y, G = K, model = "VVV",
             initialization = ini, parallel = FALSE, verbose = FALSE,
             iter.max = max(max_iter, COMPETITOR_DEFAULTS$CNmixt$iter.max)),
             error = function(e) NULL))
      if (is.null(f)) next
      keep(f, tryCatch(as.numeric(f$models[[1]]$loglik), error = function(e) -Inf))
    }
    if (is.null(best) || is.null(best$models[[1]]$group))
      return(list(fit = NULL, cluster = NULL))
    return(list(fit = best, cluster = best$models[[1]]$group))
  }

  if (grepl("MixSIT", code)) {
    mode  <- c(MixSIT_None = "none", MixSIT_U = "U", MixSIT_Full = "full")[[code]]
    inits <- rep(c("kmeans", "hclust", "random"), length.out = R)
    for (r in seq_len(R)) {
      f <- tryCatch(fitMixSIT_ECM(Y, K = K, constraint_mode = mode,
             max_iter = max_iter, tol = MIXSIT$tol,
             n_short = 1, short_iter = 1, init_methods = inits[r],
             lambda_U = MIXSIT$lambda_U, lambda_O = MIXSIT$lambda_O,
             lambda_magnitude = MIXSIT$lambda_magnitude,
             lambda_hetero = MIXSIT$lambda_hetero, verbose = FALSE),
             error = function(e) NULL)
      if (is.null(f) || is.null(f$cluster)) next
      keep(f, tryCatch(as.numeric(tail(f$loglik, 1)), error = function(e) -Inf))
    }
    return(list(fit = best, cluster = if (is.null(best)) NULL else best$cluster))
  }
  stop("unknown method code: ", code)
}

METHOD_NAMES <- c("GMM", "tMM", "SNMM", "STMM", "MSGHD", "CGMM",
                  "MixSIT_None", "MixSIT_U", "MixSIT_Full")

# Version the stored diagnostic separately from the fitted model results.  The first implementation compared latent column j with observed cell j directly.
# That comparison is not invariant to the likelihood-preserving axis permutations used by the estimator, and became visibly wrong for EE once its
# former cross-component tail-order restriction was removed.
WBAR_AXIS_MATCH_VERSION <- "global-observed-direction-v1"
SECTION4_WBAR_CELLS <- "S7_CellWise10"

# Exact square linear assignment by dynamic programming over column subsets. Here m = 10, so O(m 2^m) is negligible and avoids an extra package dependency.
# The result maps each row of the cost matrix to one distinct column.
solve_axis_assignment <- function(cost) {
  cost <- as.matrix(cost)
  if (!is.numeric(cost) || nrow(cost) != ncol(cost) || !nrow(cost) ||
      any(!is.finite(cost)))
    stop("axis-assignment cost must be a non-empty finite square matrix")
  n <- nrow(cost)
  if (n >= 30L) stop("axis-assignment dimension is too large for bit masks")
  bits <- bitwShiftL(1L, seq.int(0L, n - 1L))
  nstate <- bitwShiftL(1L, n)
  value <- rep(Inf, nstate)
  parent_mask <- parent_col <- rep(NA_integer_, nstate)
  value[1L] <- 0
  tol <- 1e-12

  for (mask in seq.int(0L, nstate - 2L)) {
    idx <- mask + 1L
    if (!is.finite(value[idx])) next
    row <- sum(bitwAnd(mask, bits) != 0L) + 1L
    for (col in seq_len(n)) if (bitwAnd(mask, bits[col]) == 0L) {
      new_mask <- bitwOr(mask, bits[col])
      candidate <- value[idx] + cost[row, col]
      new_idx <- new_mask + 1L
      if (candidate < value[new_idx] - tol) {
        value[new_idx] <- candidate
        parent_mask[new_idx] <- mask
        parent_col[new_idx] <- col
      }
    }
  }

  map <- integer(n)
  mask <- nstate - 1L
  for (row in n:1L) {
    idx <- mask + 1L
    col <- parent_col[idx]
    if (is.na(col)) stop("axis-assignment backtracking failed")
    map[row] <- col
    mask <- parent_mask[idx]
  }
  map
}

# Map fitted latent axes to the observed coordinate axes.  Column j of A = U' O' is the observed-space direction carrying latent coordinate j.  The
# absolute cosine makes the assignment sign invariant; the global assignment prevents two observed coordinates from claiming the same latent direction.
match_observed_axes <- function(component) {
  A <- t(component$U) %*% t(component$O)
  norms <- sqrt(colSums(A^2))
  if (any(!is.finite(norms)) || any(norms <= 0))
    stop("non-finite or zero observed-space direction")
  Q <- sweep(A, 2, norms, "/")
  similarity <- pmin(abs(Q), 1)
  map <- solve_axis_assignment(acos(similarity))
  list(map = map,
       abs_cosine = similarity[cbind(seq_len(nrow(Q)), map)])
}

# Posterior-average the scale weights only after matching axes within every component.  For observed coordinate l, map_k[l] is the latent axis of
# component k whose observed-space direction is globally closest to e_l. This is the appropriate diagnostic for the axis-aligned S5 intervention; it
# does not redefine the estimator or claim rotation-invariant cell weights.
mixsit_wbar <- function(Y, fit) {
  ex <- tryCatch(E_Step_cpp(Y, fit$params), error = function(e) NULL)
  if (is.null(ex)) return(NULL)
  W <- ex$W; Z <- ex$Z
  K <- dim(W)[3]
  d <- ncol(Y)
  out <- matrix(0, nrow(Y), d)
  maps <- cosines <- matrix(NA_real_, K, d)
  for (k in seq_len(K)) {
    mt <- tryCatch(match_observed_axes(fit$params$components[[k]]),
                   error = function(e) NULL)
    if (is.null(mt)) return(NULL)
    maps[k, ] <- mt$map
    cosines[k, ] <- mt$abs_cosine
    out <- out + Z[, k] * W[, mt$map, k, drop = FALSE][, , 1L]
  }
  list(weights = out, axis_map = maps, axis_abs_cosine = cosines,
       matching = WBAR_AXIS_MATCH_VERSION)
}

# Standard result provenance is not sufficient to invalidate only this derived diagnostic.  Require the Wbar schema version for contaminated Section 4 cells
# while leaving every unrelated cached fit reusable.
section4_rep_cache_valid <- function(file, cell) {
  if (!cache_valid(file)) return(FALSE)
  if (!cell %in% SECTION4_WBAR_CELLS) return(TRUE)
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  !is.null(x$wbar) &&
    identical(x$wbar$axis_match_version, WBAR_AXIS_MATCH_VERSION) &&
    all(vapply(x$wbar$wbar, function(z)
      is.list(z) && identical(z$matching, WBAR_AXIS_MATCH_VERSION) &&
        is.matrix(z$weights), logical(1)))
}

fit_all_models <- function(Y, K, max_iter, want_wbar = FALSE,
                           tag = NULL, rep = NULL) {
  N <- nrow(Y); d <- ncol(Y)
  clusters <- timing <- success <- converged <- ic <- wbar <- list()
  dropped <- list()
  for (method in METHOD_NAMES) {
    t0 <- Sys.time()
    # tag/rep give this method its own stream; see fit_multistart
      res <- tryCatch(fit_multistart(method, Y, K, R_RESTART, max_iter,
                                   tag = tag, rep = rep),
                    error = function(e) NULL)
    timing[[method]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok <- !is.null(res) && !is.null(res$cluster) && !all(is.na(res$cluster))
    dropped[[method]] <- if (ok) sum(is.na(res$cluster)) else NA_integer_
    converged[[method]] <- if (!ok) FALSE else if (method %in% c("tMM", "SNMM", "STMM")) {
      isTRUE(res$fit$error == 0)
    } else if (method == "MSGHD") {
      length(res$fit@loglik) < COMPETITOR_DEFAULTS$MSGHD$max.iter
    } else if (method == "CGMM") {
      isTRUE(res$fit$models[[1]]$iter.stop < COMPETITOR_DEFAULTS$CNmixt$iter.max)
    } else if (grepl("MixSIT", method)) {
      isTRUE(res$fit$converged)
    } else NA
    if (!ok) {
      clusters[[method]] <- rep(NA_integer_, N); success[[method]] <- FALSE
      ic[[method]] <- list(AIC = NA_real_, BIC = NA_real_, ICL = NA_real_)
    } else {
      clusters[[method]] <- res$cluster; success[[method]] <- TRUE
      ic[[method]] <- calc_information_criteria(Y, res$fit, method, K, d)
      if (want_wbar && grepl("MixSIT", method))
        wbar[[method]] <- mixsit_wbar(Y, res$fit)
    }
  }
  list(clusters = clusters, timing = timing, success = success,
       converged = converged, ic = ic, wbar = wbar, dropped = dropped)
}

# One replicate
run_one_replicate <- function(sim_id, cell, N, d, K, global_seed, max_iter) {
  rng_fix()
  data <- generate_cell(cell, N, d, K, mk_seed(global_seed, cell, "data", sim_id))
  want_wbar <- !is.null(data$contam)
  fits <- fit_all_models(data$Y, K, max_iter, want_wbar = want_wbar,
                         tag = paste(global_seed, cell, sep = "/"), rep = sim_id)
  row <- list(sim_id = sim_id)
  for (m in METHOD_NAMES) {
    row[[paste0(m, "_ARI")]]     <- calc_ari(data$labels, fits$clusters[[m]])
    row[[paste0(m, "_NMI")]]     <- calc_nmi(data$labels, fits$clusters[[m]])
    row[[paste0(m, "_AIC")]]     <- fits$ic[[m]]$AIC
    row[[paste0(m, "_BIC")]]     <- fits$ic[[m]]$BIC
    row[[paste0(m, "_ICL")]]     <- fits$ic[[m]]$ICL
    row[[paste0(m, "_Time")]]    <- fits$timing[[m]]
    row[[paste0(m, "_Success")]] <- fits$success[[m]]
    row[[paste0(m, "_Converged")]] <- fits$converged[[m]]
    row[[paste0(m, "_Dropped")]] <- fits$dropped[[m]]
  }
  out <- as.data.frame(row, stringsAsFactors = FALSE)
  # The latent scales travel as an attribute so the results table keeps its
  # one-row-per-replicate shape; the driver peels them off into its own file.
  if (want_wbar)
    attr(out, "wbar") <- list(
      sim_id = sim_id, contam = data$contam, cluster_true = data$labels,
      axis_match_version = WBAR_AXIS_MATCH_VERSION, wbar = fits$wbar)
  out
}
