# Section6_LambdaSensitivity.R
#
#   Rscript Section6_LambdaSensitivity.R <lam> <N> <rep>   one fit
#   Rscript Section6_LambdaSensitivity.R [n_parallel]      driver
#
# The supplementary experiment for Section 6.1: how sensitive is anything to the one tuning constant the model carries, lambda_het, and should it be a fixed
# number or chosen from the data?


SIM_DIR <- normalizePath(".", winslash = "/")
setwd(SIM_DIR)

LAMS  <- c(0, 0.5, 1, 2, 5, 10, 25)
NGRID <- c(200, 400, 800, 1600)
NREP  <- 10
MDIM  <- 10
KTRUE <- 3

lam_tag  <- function(l) sprintf("L%s", sub("[.]", "p", format(l, trim = TRUE)))
out_path <- function(l, N, r)
  file.path("results", sprintf("_lam_%s_%d_%02d.rds", lam_tag(l), N, r))

run_one <- function(lam, N, r) {
  suppressWarnings(suppressMessages({
    source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R") }))
  source("Section6_FlipLib.R"); source("Section6_TruthLib.R")
  source("Section6_PenaltyLib.R")
  d <- MDIM; K <- KTRUE
  install_penalty("log", lam)                 # gamma penalty off: C++ CM-2 kept
  rng_fix()
  # the SAME truth and the SAME data as Section6_Recovery.R, keyed identically,
  # so the lambda sweep is a sweep on the deployment truth rather than on a
  # parallel one that only looks like it
  TR <- make_truth("D3", d, K, mk_seed("Recovery", "D3", "truth", r))
  set.seed(mk_seed("Recovery", "D3", "data", N, r))
  n_k <- as.vector(rmultinom(1, N, rep(1/K, K)))
  Y <- matrix(0, N, d); lab <- integer(N); idx <- 1
  for (k in 1:K) if (n_k[k] > 0) {
    Y[idx:(idx + n_k[k] - 1), ] <- assemble(n_k[k], TR$mu[[k]], TR$U, TR$O[[k]],
                                            TR$nu[[k]], TR$gm[[k]])
    lab[idx:(idx + n_k[k] - 1)] <- k; idx <- idx + n_k[k] }
  set.seed(mk_seed("Recovery", "D3", "perm", N, r))
  p <- sample.int(N); Y <- Y[p, , drop = FALSE]; lab <- lab[p]

  seed_for("Lambda", "fit", lam, N, r)
  t0 <- Sys.time()
  f <- fit_with_flips(Y, K, lam, mode = "U")
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  R <- score(f$params, f$cluster, lab, TR, d, K, "D3", N, r)
  R$lam <- lam
  R$ari <- calc_ari(lab, f$cluster)
  R$nmi <- calc_nmi(lab, f$cluster)
  R$loglik <- f$loglik_obs                     # raw, comparable across lambda
  R$bic <- -2 * f$loglik_obs + count_parameters("MixSIT_U", K, d) * log(N)
  R$n_flips <- f$n_flips; R$secs <- secs
  R
}

args <- commandArgs(trailingOnly = TRUE)
source("00_RNG.R")
if (length(args) >= 3) {
  lam <- as.numeric(args[1]); N <- as.integer(args[2]); r <- as.integer(args[3])
  f <- out_path(lam, N, r)
  if (!cache_valid(f)) { z <- tryCatch(run_one(lam, N, r), error = function(e) NULL)
                         if (!is.null(z)) save_result(z, f) }
  quit(save = "no")
}

suppressWarnings(suppressMessages({ source("00_CONFIG.R"); source("Section4_Lib.R") }))
N_PAR <- if (length(args) >= 1) as.integer(args[1]) else 8
grid <- expand.grid(lam = LAMS, N = NGRID, rep = seq_len(NREP))
grid <- grid[order(grid$rep, grid$N), ]
done <- mapply(function(l, N, r) cache_valid(out_path(l, N, r)),
               grid$lam, grid$N, grid$rep)
todo <- grid[!done, , drop = FALSE]
cat(sprintf("lambda sweep: %d of %d outstanding, %d at a time\n",
            nrow(todo), nrow(grid), N_PAR))
run_process_queue(todo, "Section6_LambdaSensitivity.R",
  args_fun = function(j) c(j$lam, j$N, j$rep),
  done_fun = function(j) cache_valid(out_path(j$lam, j$N, j$rep)),
  label_fun = function(j) sprintf("lambda_%s_%d_%02d", lam_tag(j$lam), j$N, j$rep),
  n_parallel = N_PAR, timeout_s = 7200)
fs <- mapply(out_path, grid$lam, grid$N, grid$rep)
if (!all(vapply(fs, cache_valid, logical(1)))) stop("incomplete lambda grid")
D <- do.call(rbind, lapply(fs, readRDS))
save_result(D, file.path(RESULTS_DIR, "Section6_LambdaSensitivity.rds"))
cat(sprintf("collected %d of %d cells\n", length(fs), nrow(grid)))
