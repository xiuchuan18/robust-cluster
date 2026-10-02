# Section6_Recovery.R
#
#   Rscript Section6_Recovery.R <N> <variant> <rep>   one fit
#   Rscript Section6_Recovery.R [n_parallel]          driver
#
# Section 6.1: can the PARAMETERS be recovered, not just the partition?

SIM_DIR <- normalizePath(".", winslash = "/")
setwd(SIM_DIR)

NGRID    <- c(200, 400, 800, 1600, 3200)
NREP     <- 50
MDIM     <- 10
KTRUE    <- 3
DESIGN   <- "D3"
VARIANTS <- c(EV = "MixSIT_U", VV = "MixSIT_None", EE = "MixSIT_Full")

out_path <- function(N, v, r)
  file.path("results", sprintf("_rec_%s_%04d_%02d.rds", v, N, r))

run_one <- function(N, vname, r) {
  suppressWarnings(suppressMessages({
    source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R")
    source("Section6_TruthLib.R") }))
  rng_fix()
  d <- MDIM; K <- KTRUE

  # the truth is drawn from its own stream and does NOT depend on N, so the sample-size sweep is a sweep on one truth rather than on a new truth each time
  TR <- make_truth(DESIGN, d, K, mk_seed("Recovery", DESIGN, "truth", r))

  set.seed(mk_seed("Recovery", DESIGN, "data", N, r))
  n_k <- as.vector(rmultinom(1, N, rep(1 / K, K)))
  Y <- matrix(0, N, d); lab <- integer(N); i <- 1
  for (k in seq_len(K)) if (n_k[k] > 0) {
    Y[i:(i + n_k[k] - 1), ] <- assemble(n_k[k], TR$mu[[k]], TR$U, TR$O[[k]],
                                        TR$nu[[k]], TR$gm[[k]])
    lab[i:(i + n_k[k] - 1)] <- k; i <- i + n_k[k] }
  set.seed(mk_seed("Recovery", DESIGN, "perm", N, r))
  p <- sample.int(N); Y <- Y[p, , drop = FALSE]; lab <- lab[p]

  code <- VARIANTS[[vname]]
  t0 <- Sys.time()
  fit <- tryCatch(fit_multistart(code, Y, K, R_RESTART, MIXSIT$max_iter,
                                 tag = paste("Recovery", DESIGN, sep = "/"),
                                 rep = paste(N, r, sep = "_")),
                  error = function(e) NULL)
  if (is.null(fit) || is.null(fit$cluster)) return(NULL)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  R <- tryCatch(score(fit$fit$params, fit$cluster, lab, TR, d, K, DESIGN, N, r),
                error = function(e) NULL)
  if (is.null(R)) return(NULL)
  # VV carries a U per component; score() reads component 1 only, which is right for the shared-U variants and wrong here
  if (code == "MixSIT_None")
    R$U_err <- mean(vapply(fit$fit$params$components,
                           function(cp) sqrt(sum((cp$U - TR$U)^2)) / d, 0))
  R$variant <- vname
  R$ari     <- calc_ari(lab, fit$cluster)
  R$secs    <- secs
  R$iters   <- fit$fit$iterations
  R$n_flips <- if (is.null(fit$fit$n_flips)) NA_integer_ else fit$fit$n_flips
  R$loglik  <- fit$fit$loglik_obs
  R
}

args <- commandArgs(trailingOnly = TRUE)
source("00_RNG.R")
if (length(args) >= 3) {
  N <- as.integer(args[1]); v <- args[2]; r <- as.integer(args[3])
  f <- out_path(N, v, r)
  if (!cache_valid(f)) {
    z <- tryCatch(run_one(N, v, r), error = function(e) NULL)
    if (!is.null(z)) save_result(z, f)
  }
  quit(save = "no")
}

suppressWarnings(suppressMessages({ source("00_CONFIG.R"); source("Section4_Lib.R") }))
N_PAR <- if (length(args) >= 1) as.integer(args[1]) else 8
grid <- expand.grid(N = NGRID, variant = names(VARIANTS), rep = seq_len(NREP),
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$rep, grid$N), ]
done <- mapply(function(N, v, r) cache_valid(out_path(N, v, r)),
               grid$N, grid$variant, grid$rep)
todo <- grid[!done, , drop = FALSE]
cat(sprintf("recovery study: %d of %d fits outstanding, %d at a time\n",
            nrow(todo), nrow(grid), N_PAR))
run_process_queue(todo, "Section6_Recovery.R",
  args_fun = function(j) c(j$N, j$variant, j$rep),
  done_fun = function(j) cache_valid(out_path(j$N, j$variant, j$rep)),
  label_fun = function(j) sprintf("recovery_%s_%04d_%02d", j$variant, j$N, j$rep),
  n_parallel = N_PAR, timeout_s = 7200)
fs <- mapply(out_path, grid$N, grid$variant, grid$rep)
if (!all(vapply(fs, cache_valid, logical(1)))) stop("incomplete recovery grid")
D <- do.call(rbind, lapply(fs, readRDS))
save_result(D, file.path(RESULTS_DIR, "Section6_Recovery.rds"))
cat(sprintf("collected %d of %d fits\n", length(fs), nrow(grid)))
