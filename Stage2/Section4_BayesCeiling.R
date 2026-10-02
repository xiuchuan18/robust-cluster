# Section4_BayesCeiling.R
# Output: results/Section4_BayesCeiling.{rds,csv} and a printed table.

source("00_CONFIG.R")
source("Section4_Lib.R")          # generators, verbatim, no simulation is run
suppressPackageStartupMessages({
  library(mvtnorm); library(mclust); library(EMMIXuskew)
})


# Component log-densities, one per generator

# S1, S2: plain Gaussian components
logdens_gauss <- function(Y, mu_list, Sigma)
  sapply(mu_list, function(mu) mvtnorm::dmvnorm(Y, mu, Sigma, log = TRUE))

# S5: N(mu, I) components with additive per-cell contamination, so each
# cell is (1-alpha) N(mu, 1) + alpha N(mu, 1 + eta^2), independent across cells.
logdens_cellcontam <- function(Y, mu_list, alpha, eta) {
  d <- ncol(Y); s2 <- sqrt(1 + eta^2)
  sapply(mu_list, function(mu) {
    lp <- numeric(nrow(Y))
    for (j in seq_len(d))
      lp <- lp + log((1 - alpha) * dnorm(Y[, j], mu[j], 1) +
                     alpha       * dnorm(Y[, j], mu[j], s2))
    lp
  })
}

# S3: Gaussian copula with t margins, closed form:
#   x_j = s_j (y_j - mu_j), s_j = sqrt(nu_j / (nu_j - 2));  u_j = F_{t_nu_j}(x_j);
#   z_j = Phi^-1(u_j);  f(y) = [phi_R(z) / prod_j phi(z_j)] prod_j f_{t_nu_j}(x_j) s_j
logdens_gc <- function(Y, mu_list, R, nuj_list) {
  sapply(seq_along(mu_list), function(k) {
    Y <- as.matrix(Y); d <- ncol(Y); nuj <- nuj_list[[k]][1:d]
    sg <- sqrt(nuj / (nuj - 2))
    X  <- sweep(sweep(Y, 2, mu_list[[k]]), 2, sg, "*")
    U  <- vapply(seq_len(d), function(j) stats::pt(X[, j], df = nuj[j]), numeric(nrow(Y)))
    U  <- pmin(pmax(U, 1e-14), 1 - 1e-14)
    Z  <- stats::qnorm(U)
    lc <- mvtnorm::dmvnorm(Z, mean = rep(0, d), sigma = R, log = TRUE) -
          rowSums(stats::dnorm(Z, log = TRUE))
    lm <- rowSums(vapply(seq_len(d), function(j)
            stats::dt(X[, j], df = nuj[j], log = TRUE), numeric(nrow(Y))))
    lc + lm + sum(log(sg))
  })
}

# S3, S4: unrestricted multivariate skew-t, the family EMMIXuskew::rmst samples
logdens_umst <- function(Y, mu_list, sigma_list, delta_list, dof)
  sapply(seq_along(mu_list), function(k)
    log(pmax(as.numeric(EMMIXuskew::dmst(Y, mu_list[[k]], sigma_list[[k]],
                                         delta_list[[k]], dof)), 1e-300)))

# Fernandez-Steel skew-t log-density, the source rfs() actually samples:
#   f(e) = 2/(g + 1/g) * t_nu(e/g)   e >= 0
#   f(e) = 2/(g + 1/g) * t_nu(e*g)   e <  0
logdens_fs <- function(e, nu, g)
  log(2) - log(g + 1 / g) +
  ifelse(e >= 0, dt(e / g, nu, log = TRUE), dt(e * g, nu, log = TRUE))

# S5: y = mode_k + O_k' eps with eps_j independent FS(nu_j, gamma_j) and
# mode_k = mu_k - O_k' E[eps]. |det O| = 1, so there is no Jacobian term.
fs_mode <- function(mu, O, nu, gm) {
  me <- vapply(seq_along(nu), function(j) fs_mean(gm[j], nu[j]), numeric(1))
  mu - as.vector(t(O) %*% me)
}
logdens_S5 <- function(Y, mu_list, O_list, nu_list, gm_list) {
  d <- ncol(Y)
  sapply(seq_along(mu_list), function(k) {
    mode_k <- fs_mode(mu_list[[k]], O_list[[k]], nu_list[[k]], gm_list[[k]])
    E <- sweep(Y, 2, mode_k, "-") %*% t(O_list[[k]])       # eps = O (y - mode)
    lp <- numeric(nrow(Y))
    for (j in seq_len(d))
      lp <- lp + logdens_fs(E[, j], nu_list[[k]][j], gm_list[[k]][j])
    lp
  })
}


truth_for <- function(cell, N, d, seed) {
  p3 <- rep(1/3, 3)

  if (cell == "S1_Gaussian") {
    mus <- list(rep(0, d), c(4, 4, rep(0, d - 2)), c(4, -2, rep(0, d - 2)))
    return(list(pi = p3, f = function(Y) logdens_gauss(Y, mus, diag(d))))
  }
  if (cell == "S2_Imbalanced") {
    mus <- list(rep(0, d), c(2, 2, rep(0, d - 2)), c(-2, -2, rep(0, d - 2)))
    return(list(pi = c(0.1, 0.3, 0.6),
                f = function(Y) logdens_gauss(Y, mus, diag(d))))
  }
  if (cell == "S7_CellWise10") {
    h <- max(2, min(4, floor(d / 2)))
    mus <- list(numeric(d), mu_at(d, 1:h, 2.5), mu_at(d, (h + 1):(2 * h), 2.5))
    return(list(pi = p3,
                f = function(Y) logdens_cellcontam(Y, mus, 0.05, 10)))
  }
  if (cell == "S3_GaussCopula") {
    mus <- lapply(1:3, function(k) S3_SCALE * mu_at(d, 1:2, S3_CENTRE[[k]]))
    R <- matrix(S3_GC_RHO, d, d); diag(R) <- 1
    return(list(pi = p3, f = function(Y) logdens_gc(Y, mus, R, S3_NU_AXIS)))
  }
  if (cell == "S6_AlignedWide") {
    set.seed(seed); invisible(rmultinom(1, N, p3))
    Os  <- lapply(1:3, function(k) qr.Q(qr(matrix(rnorm(d^2), d, d))))
    mus <- list(rep(0, d), mu_at(d, 1:3, 2), mu_at(d, 4:6, 2))
    return(list(pi = p3, f = function(Y) logdens_S5(Y, mus, Os, S6_NU(d), S6_GM(d))))
  }
  if (cell == "S4_Skewed") {
    mus <- list(rep(0, d), c(2, 2, rep(0, d - 2)), c(-2, -2, rep(0, d - 2)))
    del <- list(c(10,6,5,-8,4)[1:d], c(12,4,8,-6,6)[1:d], c(10,-5,-4,-4,10)[1:d])
    sg  <- list(get_fixed_covariance(d, 1), get_fixed_covariance(d, 2),
                get_fixed_covariance(d, 3))
    return(list(pi = p3, f = function(Y) logdens_umst(Y, mus, sg, del, 20)))
  }
  stop("no truth defined for cell: ", cell)
}


N_REP      <- if (exists("BAYES_N_REP")) BAYES_N_REP else 100

rows <- list()
for (cell in names(CONFIG_SIM$scenarios)) {
  sc <- CONFIG_SIM$scenarios[[cell]]
  N <- sc$N; d <- sc$d
  mc <- FALSE
  nr <- N_REP
  ari <- conf <- amb <- numeric(nr)
  t0 <- Sys.time()
  for (r in seq_len(nr)) {
    # the SAME data seed run_one_replicate() uses (00_RNG.R role streams), so
    # the ceiling is computed on exactly the 100 samples the methods were fitted to
    seed <- mk_seed(CONFIG_SIM$seed, cell, "data", r)
    dat  <- generate_cell(cell, N, d, 3, seed)
    tr   <- truth_for(cell, N, d, seed)
    L    <- sweep(tr$f(dat$Y), 2, log(tr$pi), "+")
    mx   <- apply(L, 1, max)
    Z    <- exp(L - mx); Z <- Z / rowSums(Z)
    ari[r]  <- mclust::adjustedRandIndex(dat$labels, max.col(Z, "first"))
    conf[r] <- mean(apply(Z, 1, max))
    amb[r]  <- mean(apply(Z, 1, max) < 0.9)
  }
  rows[[cell]] <- data.frame(
    cell = cell, N = N, d = d, n_rep = nr, monte_carlo = mc,
    ARI_Bayes = mean(ari), ARI_Bayes_sd = sd(ari),
    ARI_Bayes_se = sd(ari) / sqrt(nr),
    mean_max_post = mean(conf), pct_ambiguous = 100 * mean(amb))
  cat(sprintf("%-16s ARI_Bayes = %.4f (sd %.4f, se %.4f)%s  max-post %.3f  ambiguous %.1f%%  [%.1f min]\n",
              cell, mean(ari), sd(ari), sd(ari) / sqrt(nr),
              if (mc) " [MC]" else "    ", mean(conf), 100 * mean(amb),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
ceil <- do.call(rbind, rows)
save_result(ceil, file.path(RESULTS_DIR, "Section4_BayesCeiling.rds"))
write_csv_atomic(ceil, file.path(RESULTS_DIR, "Section4_BayesCeiling.csv"))

# attainment table, if the simulation results are on disk
res_file <- file.path(RESULTS_DIR, "Section4_Simulation_Results.rds")
if (file.exists(res_file)) {
  res <- readRDS(res_file)
  cat("\n===== attained fraction of the ceiling =====\n")
  cat(sprintf("%-16s %10s %10s %12s %12s\n",
              "cell", "ceiling", "best ARI", "best method", "attained"))
  for (cell in names(CONFIG_SIM$scenarios)) {
    if (is.null(res[[cell]])) next
    cols  <- grep("_ARI$", names(res[[cell]]), value = TRUE)
    means <- vapply(cols, function(cc) mean(res[[cell]][[cc]], na.rm = TRUE), 0)
    b   <- which.max(means)
    lab <- METHOD_CODE_TO_LABEL[sub("_ARI$", "", names(means)[b])]
    cb  <- ceil$ARI_Bayes[ceil$cell == cell]
    cat(sprintf("%-16s %10.4f %10.4f %12s %11.1f%%\n",
                cell, cb, means[b], lab, 100 * means[b] / cb))
  }
  cat("\nA method at 90-100% of the ceiling is not underperforming; the cell\n")
  cat("simply does not contain enough separation for a higher ARI to exist.\n")
} else {
  cat("\n(Section4_Simulation_Results.rds not found; attainment table skipped.)\n")
}
