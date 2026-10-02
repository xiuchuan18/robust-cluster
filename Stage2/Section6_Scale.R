# Section6_Scale.R
#
#   Rscript Section6_Scale.R <N> <m> <K> <method> <rep>   one timed fit
#   Rscript Section6_Scale.R [n_parallel]                 driver
#
# Section 6.3: how the cost of a MixSIT fit moves with sample size, dimension and
# the number of components, against a cheap baseline and the flexible competitor.
#
# Three one-factor sweeps through a common anchor at (N = 800, m = 10, K = 3), so
# each exponent is read at a fixed level of the other two and the anchor is
# measured once and shared:
#     N in {200, 400, 800, 1600}   at m = 10, K = 3
#     m in {3, 5, 10, 15, 20}      at N = 800, K = 3
#     K in {2, 3, 4, 5, 6}         at N = 800, m = 10

SIM_DIR <- normalizePath(".", winslash = "/")
setwd(SIM_DIR)

NREP    <- 30   # 10 -> 30 on 2026-09-14; replicates 1-10 are the cached ones
METHODS <- c("GMM", "MSGHD", "MixSIT_None", "MixSIT_U", "MixSIT_Full")
ANCHOR  <- list(N = 800, m = 10, K = 3)
SWEEPS  <- list(
  N = lapply(c(200, 400, 800, 1600), function(x) list(N = x, m = 10, K = 3)),
  m = lapply(c(3, 5, 10, 15, 20),    function(x) list(N = 800, m = x, K = 3)),
  K = lapply(c(2, 3, 4, 5, 6),       function(x) list(N = 800, m = 10, K = x)))

out_path <- function(N, m, K, meth, r)
  file.path("results", sprintf("_scale_N%d_m%d_K%d_%s_%02d.rds", N, m, K, meth, r))

run_one <- function(N, m, K, method, r) {
  suppressWarnings(suppressMessages({
    source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R") }))
  source("Section6_FlipLib.R"); source("Section6_ScaleLib.R")
  suppressPackageStartupMessages(library(MixGHD))
  rng_fix()
  d <- gen_scale(N, m, K, mk_seed("Scale", "data", N, m, K, r))
  z <- time_fit(method, d$Y, K, MIXSIT$lambda_hetero,
                tag = paste("Scale", N, m, K, sep = "/"), rep = r)
  data.frame(N = N, m = m, K = K, method = method, rep = r,
             t_base = z$t_base, t_flip = z$t_flip, t_tot = z$t_base + z$t_flip,
             iters = z$iters, flips = z$flips, ok = z$ok,
             stringsAsFactors = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
source("00_RNG.R")
if (length(args) >= 5) {
  N <- as.integer(args[1]); m <- as.integer(args[2]); K <- as.integer(args[3])
  f <- out_path(N, m, K, args[4], as.integer(args[5]))
  if (!cache_valid(f)) {
    z <- tryCatch(run_one(N, m, K, args[4], as.integer(args[5])),
                  error = function(e) NULL)
    if (!is.null(z)) save_result(z, f) }
  quit(save = "no")
}

suppressWarnings(suppressMessages({ source("00_CONFIG.R"); source("Section4_Lib.R") }))
N_PAR <- if (length(args) >= 1) as.integer(args[1]) else 8

# distinct configurations, so the anchor is not fitted three times
cfg <- unique(do.call(rbind, lapply(unlist(SWEEPS, recursive = FALSE),
  function(s) data.frame(N = s$N, m = s$m, K = s$K))))
grid <- do.call(rbind, lapply(seq_len(NREP), function(r)
  do.call(rbind, lapply(METHODS, function(mm)
    data.frame(cfg, method = mm, rep = r, stringsAsFactors = FALSE)))))
# cheap methods first inside a replicate so a partial read is usable
grid$prio <- match(grid$method, c("GMM", "MixSIT_Full", "MixSIT_U",
                                  "MixSIT_None", "MSGHD"))
grid <- grid[order(grid$rep, grid$prio, grid$m, grid$N, grid$K), ]
done <- mapply(function(N, m, K, meth, r)
  cache_valid(out_path(N, m, K, meth, r)), grid$N, grid$m, grid$K,
  grid$method, grid$rep)
todo <- grid[!done, , drop = FALSE]
cat(sprintf("scalability grid: %d of %d timed fits outstanding, %d at a time\n",
            nrow(todo), nrow(grid), N_PAR))
run_process_queue(todo, "Section6_Scale.R",
  args_fun = function(j) c(j$N, j$m, j$K, j$method, j$rep),
  done_fun = function(j) cache_valid(out_path(j$N, j$m, j$K, j$method, j$rep)),
  label_fun = function(j) sprintf("scale_%d_%d_%d_%s_%02d", j$N, j$m,
                                  j$K, j$method, j$rep),
  n_parallel = N_PAR, timeout_s = 7200)
fs <- mapply(out_path, grid$N, grid$m, grid$K, grid$method, grid$rep)
if (!all(vapply(fs, cache_valid, logical(1)))) stop("incomplete scaling grid")
D <- do.call(rbind, lapply(fs, readRDS))
save_result(D, file.path(RESULTS_DIR, "Section6_Scale.rds"))
cat(sprintf("collected %d of %d timed fits\n", length(fs), nrow(grid)))
