# Confirmatory phase map: eight frozen cells, 50 fresh replicates, K=1,...,5.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
SCRIPT <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
STAGE2 <- dirname(SCRIPT)
setwd(STAGE2)
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R")
  source("Section5_Datasets.R"); source("CalFlex_FormalLib.R")
}))

OUT <- file.path(RESULTS_DIR, "calibrated-flexibility")
CACHE <- file.path(OUT, "cache", "phase-map")
LOG <- file.path(OUT, "logs", "phase-map")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

KMAX <- 5L; N <- 400L; D <- 5L; KTRUE <- 3L; NREP <- 50L
CELLS <- data.frame(
  mechanism = c("cellSD", "cellSD", "cellSD", "cellSD", "cellSD", "cellSD",
                "tailMM", "skewMM"),
  severity = c(0, 2/3, 1, 0, 2/3, 1, 1, 1),
  separation = c(4, 4, 4, 2.5, 2.5, 2.5, 2.5, 2.5),
  stringsAsFactors = FALSE
)

harmonic_path <- function(base, target, s) 1 / ((1 - s) / base + s / target)
balanced_counts <- function(n, k) {
  z <- rep(n %/% k, k)
  if (sum(z) < n) z[seq_len(n - sum(z))] <- z[seq_len(n - sum(z))] + 1L
  z
}
phase_centres <- function(d, sep) list(
  rep(0, d), c(sep, sep, rep(0, d - 2)),
  c(sep, -sep, rep(0, d - 2)))
moment_match_identity <- function(x) {
  xc <- sweep(x, 2, colMeans(x), "-")
  S <- cov(xc) + diag(1e-8, ncol(xc))
  xc %*% solve(chol(S))
}

gen_phase_formal <- function(mechanism, severity, sep, rep,
                             n = N, d = D, k = KTRUE) {
  seed <- calflex_seed("PhaseMap", mechanism, severity, sep, "data", rep)
  set.seed(seed)
  nk <- balanced_counts(n, k); mus <- phase_centres(d, sep)
  y <- matrix(NA_real_, n, d); lab <- rep(seq_len(k), nk); at <- 1L
  if (mechanism == "tailMM") {
    targets <- list(c(3,5,8,15,30), c(30,15,8,5,3), c(8,30,3,15,5))
    R <- matrix(0.30, d, d); diag(R) <- 1; L <- chol(R)
    for (j in seq_len(k)) {
      nn <- nk[j]; nu <- harmonic_path(200, targets[[j]][seq_len(d)], severity)
      u <- pnorm(matrix(rnorm(nn * d), nn, d) %*% L)
      u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
      x <- vapply(seq_len(d), function(q)
        qt(u[, q], df = nu[q]) / sqrt(nu[q] / (nu[q] - 2)), numeric(nn))
      x <- moment_match_identity(x)
      rows <- at:(at + nn - 1L); y[rows, ] <- sweep(x, 2, mus[[j]], "+")
      at <- at + nn
    }
  } else if (mechanism == "skewMM") {
    targets <- list(c(5,3,2,-3,2), c(5,2,3,-2,3), c(4,-2,-2,-2,4))
    nu <- harmonic_path(200, 20, severity)
    for (j in seq_len(k)) {
      nn <- nk[j]; delta <- severity * targets[[j]][seq_len(d)]
      x <- rmst_sdb(nn, rep(0, d), diag(d), delta, nu)
      rows <- at:(at + nn - 1L)
      y[rows, ] <- sweep(moment_match_identity(x), 2, mus[[j]], "+")
      at <- at + nn
    }
  } else if (mechanism == "cellSD") {
    for (j in seq_len(k)) {
      nn <- nk[j]; rows <- at:(at + nn - 1L)
      x <- mvtnorm::rmvnorm(nn, rep(0, d), diag(d))
      mask <- matrix(runif(nn * d) < 0.05, nn, d)
      shock <- matrix(rnorm(nn * d), nn, d)
      x[mask] <- x[mask] + severity * 8 * shock[mask]
      x <- scale(x)
      y[rows, ] <- sweep(x, 2, mus[[j]], "+"); at <- at + nn
    }
  } else stop("unknown frozen phase mechanism: ", mechanism)
  p <- sample.int(n)
  list(Y = scale(y[p, , drop = FALSE]), labels = lab[p], data_seed = seed)
}

phase_cache <- function(mechanism, severity, sep, method, rep)
  file.path(CACHE, sprintf("phase_%s_s%03d_d%03d_%s_r%02d.rds", mechanism,
    round(100 * severity), round(10 * sep), method, rep))

fit_phase_path <- function(mechanism, severity, sep, method, rep) {
  dat <- gen_phase_formal(mechanism, severity, sep, rep)
  rows <- vector("list", KMAX)
  for (kk in seq_len(KMAX)) {
    fit_tag <- paste(CALFLEX_SEED_PREFIX, "PhaseMap", mechanism, severity,
                     sep, sep = "/")
    fit_rep <- paste(method, rep, "K", kk, sep = "_")
    fit_seed <- calflex_fit_seed(fit_tag, method, fit_rep)
    set.seed(fit_seed)
    t0 <- proc.time()[3]
    rr <- tryCatch(fit_multistart(method, dat$Y, kk, R_RESTART,
      COMPETITOR_MAX_ITER, tag = fit_tag, rep = fit_rep), error = function(e) NULL)
    elapsed <- proc.time()[3] - t0
    ok <- !is.null(rr) && !is.null(rr$fit) && !is.null(rr$cluster)
    ll <- if (ok) calc_loglik(dat$Y, rr$fit, method) else NA_real_
    p_actual <- if (ok) calflex_npar(method, kk, D) else NA_real_
    p_package <- if (method == "MSGHD") calflex_msghd_npar(kk, D, "package") else p_actual
    p_over <- if (method == "MSGHD") calflex_msghd_npar(kk, D, "over") else p_actual
    ent <- if (ok) calflex_entropy(calflex_posterior(rr$fit, method)) else NA_real_
    rows[[kk]] <- data.frame(
      mechanism, severity, separation = sep, method, rep, K_true = KTRUE,
      K = kk, n = N, d = D, ok = ok && is.finite(ll),
      converged = calflex_converged(rr, method), loglik = ll,
      npar = p_actual, npar_package = p_package, npar_over = p_over,
      BIC = -2 * ll + p_actual * log(N),
      BIC_package = -2 * ll + p_package * log(N),
      BIC_over = -2 * ll + p_over * log(N),
      ICL = -2 * ll + p_actual * log(N) + 2 * ent,
      ARI = if (ok) calc_ari(dat$labels, rr$cluster) else NA_real_,
      seconds = elapsed, data_seed = dat$data_seed, fit_seed = fit_seed)
  }
  do.call(rbind, rows)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 5L && args[1] != "driver") {
  mechanism <- args[1]; severity <- as.numeric(args[2]); sep <- as.numeric(args[3])
  method <- args[4]; rep <- as.integer(args[5])
  f <- phase_cache(mechanism, severity, sep, method, rep)
  if (!calflex_cache_valid(f))
    calflex_save(fit_phase_path(mechanism, severity, sep, method, rep), f)
  quit(save = "no")
}
if (!length(args) || args[1] != "driver")
  stop("use five job arguments or: driver [parallel]")

np <- if (length(args) >= 2L) as.integer(args[2]) else 8L
grid <- merge(CELLS, expand.grid(method = calflex_phase_methods,
  rep = seq_len(NREP), stringsAsFactors = FALSE))
files <- mapply(phase_cache, grid$mechanism, grid$severity, grid$separation,
                grid$method, grid$rep, SIMPLIFY = FALSE)
todo <- grid[!vapply(files, calflex_cache_valid, logical(1)), , drop = FALSE]
cat(sprintf("formal phase map: %d/%d paths outstanding; %d concurrent\n",
            nrow(todo), nrow(grid), np)); flush.console()
if (nrow(todo)) {
  active <- list(); next_job <- 1L
  launch <- function(i) {
    j <- todo[i, ]; log <- file.path(LOG, sprintf(
      "%s_s%03d_d%03d_%s_r%02d.log", j$mechanism, round(100*j$severity),
      round(10*j$separation), j$method, j$rep))
    processx::process$new(file.path(R.home("bin"), "Rscript.exe"),
      c(SCRIPT, j$mechanism, j$severity, j$separation, j$method, j$rep),
      stdout = log, stderr = "2>&1", cleanup = TRUE,
      env = c(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"))
  }
  while (next_job <= nrow(todo) || length(active)) {
    while (next_job <= nrow(todo) && length(active) < np) {
      active[[as.character(next_job)]] <- launch(next_job); next_job <- next_job + 1L
    }
    Sys.sleep(0.5)
    fin <- names(active)[!vapply(active, function(p) p$is_alive(), logical(1))]
    for (id in fin) {
      p <- active[[id]]
      if (p$get_exit_status() != 0L) warning("phase job failed: ", id)
      active[[id]] <- NULL
    }
  }
}
files <- mapply(phase_cache, grid$mechanism, grid$severity, grid$separation,
                grid$method, grid$rep, SIMPLIFY = FALSE)
if (!all(vapply(files, calflex_cache_valid, logical(1))))
  stop("formal phase map incomplete")
Dall <- do.call(rbind, lapply(files, readRDS))
calflex_save(Dall, file.path(OUT, "CalFlex_PhaseMap.rds"))
calflex_write_csv(Dall, file.path(OUT, "CalFlex_PhaseMap.csv"))
cat("formal phase map complete: ", nrow(Dall), " K fits\n", sep = "")
