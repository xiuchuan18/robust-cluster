# Section6_KSelect.R
#
#   Rscript Section6_KSelect.R <cell> <method> <rep>   one scenario x method
#   Rscript Section6_KSelect.R <cell> <method> <rep> <K>  one isolated fit
#   Rscript Section6_KSelect.R [n_parallel]            driver
#
# Section 6.2: choosing the number of components by AIC, BIC and ICL.

SIM_DIR <- normalizePath(".", winslash = "/")
setwd(SIM_DIR)

KMAX    <- 6
NREP    <- 30   # 10 -> 30 on 2026-09-14; replicates 1-10 are the cached ones
NOBS    <- 600
METHODS <- c("GMM", "MSGHD", "MixSIT_None", "MixSIT_U", "MixSIT_Full")

out_path <- function(cell, meth, r)
  file.path("results", sprintf("_ksel_%s_%s_%02d.rds", cell, meth, r))

part_path <- function(cell, meth, r, K)
  file.path("results", "_ksel_parts",
            sprintf("_kselpart_%s_%s_%02d_K%02d.rds", cell, meth, r, K))

source("Section6_KSelect_fitlib.R")   # fit_at_K()

run_one_k <- function(cell, method, r, K) {
  suppressWarnings(suppressMessages({
    source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R") }))
  source("Section6_KSelectLib.R"); source("Section6_KSelect_fitlib.R")
  suppressPackageStartupMessages(library(MixGHD))
  RNGkind("Mersenne-Twister")

  # role streams (00_RNG.R): the data and the row permutation are drawn from the
  # cell's own stream, and each method is seeded inside fit_multistart, so the
  # arms are paired within a replicate and no method can disturb another.
  rng_fix()
  dat <- gen_kcell(cell, NOBS, mk_seed("KSelect", cell, "data", r))
  Y <- dat$Y; d <- ncol(Y); n <- nrow(Y)
  seed_for("KSelect", cell, "perm", r)
  p <- sample.int(n); Y <- Y[p, , drop = FALSE]; lab <- dat$labels[p]

  t0 <- Sys.time()
  f <- fit_at_K(method, Y, K, d, tag = paste("KSelect", cell, sep = "/"), rep = r)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !is.null(f) && is.finite(f$ll)
  ic <- if (ok) ic_of(f$ll, f$z, method, K, n, d)
        else list(npar = NA_integer_, AIC = NA_real_, BIC = NA_real_, ICL = NA_real_)
  data.frame(
    cell = cell, method = method, rep = r, K_true = dat$K, K = K,
    n = n, m = d, loglik = if (ok) f$ll else NA_real_,
    npar = ic$npar, AIC = ic$AIC, BIC = ic$BIC, ICL = ic$ICL,
    ari = if (ok && !is.null(f$cl)) calc_ari(lab, f$cl) else NA_real_,
    secs = secs, ok = ok, stringsAsFactors = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
source("00_RNG.R")
if (length(args) >= 4) {
  f <- part_path(args[1], args[2], as.integer(args[3]), as.integer(args[4]))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  if (!cache_valid(f)) {
    z <- tryCatch(run_one_k(args[1], args[2], as.integer(args[3]), as.integer(args[4])),
                  error = function(e) {
                    message("isolated fit failed: ", conditionMessage(e)); NULL
                  })
    if (!is.null(z)) save_result(z, f)
  }
  quit(save = "no")
}

if (length(args) >= 3) {
  f <- out_path(args[1], args[2], as.integer(args[3]))
  suppressWarnings(suppressMessages({ source("00_CONFIG.R") }))
  if (!cache_valid(f)) {
    timeout_secs <- suppressWarnings(as.numeric(Sys.getenv(
      "KSELECT_FIT_TIMEOUT_SECS", unset = "1800")))
    if (!is.finite(timeout_secs) || timeout_secs <= 0) timeout_secs <- 1800
    rows <- vector("list", KMAX)
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Section6_KSelect.R"
    script <- normalizePath(script, winslash = "/", mustWork = TRUE)
    for (K in seq_len(KMAX)) {
      pf <- part_path(args[1], args[2], as.integer(args[3]), K)
      dir.create(dirname(pf), recursive = TRUE, showWarnings = FALSE)
      if (!cache_valid(pf)) {
        status <- tryCatch({
          processx::run("Rscript",
                        c(script, args[1], args[2], args[3], as.character(K)),
                        timeout = timeout_secs * 1000,
                        error_on_status = FALSE, echo = FALSE)
          "returned"
        }, error = function(e) {
          message(sprintf("isolated fit %s/%s/rep%s/K%s stopped: %s",
                          args[1], args[2], args[3], K, conditionMessage(e)))
          "failed_or_timed_out"
        })
      }
      if (cache_valid(pf)) {
        rows[[K]] <- readRDS(pf)
      } else {
        suppressWarnings(suppressMessages({
          source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R") }))
        source("Section6_KSelectLib.R")
        dat <- gen_kcell(args[1], NOBS,
                         mk_seed("KSelect", args[1], "data", as.integer(args[3])))
        rows[[K]] <- data.frame(
          cell = args[1], method = args[2], rep = as.integer(args[3]),
          K_true = dat$K, K = K, n = nrow(dat$Y), m = ncol(dat$Y),
          loglik = NA_real_, npar = NA_integer_, AIC = NA_real_,
          BIC = NA_real_, ICL = NA_real_, ari = NA_real_,
          secs = timeout_secs, ok = FALSE, stringsAsFactors = FALSE)
      }
    }
    save_result(do.call(rbind, rows), f)
  }
  quit(save = "no")
}

suppressWarnings(suppressMessages({ source("00_CONFIG.R"); source("Section4_Lib.R") }))
source("Section6_KSelectLib.R")
N_PAR <- if (length(args) >= 1) as.integer(args[1]) else 8
grid <- expand.grid(cell = KCELLS, method = METHODS, rep = seq_len(NREP),
                    stringsAsFactors = FALSE)
# cheap methods first within a replicate, so a partial read is usable early
grid$prio <- match(grid$method, c("GMM", "MixSIT_Full", "MixSIT_U",
                                  "MixSIT_None", "MSGHD"))
grid <- grid[order(grid$rep, grid$prio, grid$cell), ]
done <- mapply(function(cell, meth, r) cache_valid(out_path(cell, meth, r)),
               grid$cell, grid$method, grid$rep)
todo <- grid[!done, , drop = FALSE]
cat(sprintf("K-selection grid: %d of %d outstanding, %d at a time\n",
            nrow(todo), nrow(grid), N_PAR))
run_process_queue(todo, "Section6_KSelect.R",
  args_fun = function(j) c(j$cell, j$method, j$rep),
  done_fun = function(j) cache_valid(out_path(j$cell, j$method, j$rep)),
  label_fun = function(j) sprintf("kselect_%s_%s_%02d", j$cell, j$method, j$rep),
  n_parallel = N_PAR, timeout_s = 14400)
fs <- mapply(out_path, grid$cell, grid$method, grid$rep)
if (!all(vapply(fs, cache_valid, logical(1)))) stop("incomplete K-selection grid")
D <- do.call(rbind, lapply(fs, readRDS))
save_result(D, file.path(RESULTS_DIR, "Section6_KSelect.rds"))
cat(sprintf("collected %d of %d jobs\n", length(fs), nrow(grid)))
