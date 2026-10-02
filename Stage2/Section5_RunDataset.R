# Section5_RunDataset.R
#
#   Rscript Section5_RunDataset.R <dataset-name>
#
# Fits every method to one clean real dataset under exactly the Section 4 protocol: one complete fit per method (R_RESTART = 1); one native k-means
# initialization call with nstart = 5 for every non-GMM method; mclust's hierarchical default for GMM; method-specific iteration caps and tolerances;
# K fixed at the known number of classes; all variables standardized. Writes results/_real_<dataset>.rds.
#
# One dataset per operating-system process, for the same reason Section 4 ends up doing this: a fitting routine that kills its worker then costs only its own
# dataset. Section5_RealData.R is the driver.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript Section5_RunDataset.R <dataset>")
nm <- args[1]

setwd(normalizePath(".", winslash = "/"))
out <- file.path("results", paste0("_real_", nm, ".rds"))

suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R")
  source("Section5_Datasets.R")
}))

d <- load_dataset(nm)
cat(sprintf("%s: n=%d m=%d K=%d\n", nm, d$n, d$m, d$K)); flush.console()

rows <- list()
rng_fix()
for (meth in METHOD_NAMES) {
  # each method on its own stream (00_RNG.R); the shared set.seed() this replaces made every method's start depend on how many draws the preceding
  # methods had consumed
  t0 <- Sys.time()
  r <- tryCatch(fit_multistart(meth, d$X, d$K, R_RESTART, COMPETITOR_MAX_ITER,
                               tag = paste("Section5", nm, sep = "/"), rep = 1),
                error = function(e) NULL)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !is.null(r) && !is.null(r$cluster) && !all(is.na(r$cluster))
  conv <- if (!ok) FALSE else if (meth %in% c("tMM", "SNMM", "STMM")) {
    isTRUE(r$fit$error == 0)
  } else if (meth == "MSGHD") {
    length(r$fit@loglik) < COMPETITOR_DEFAULTS$MSGHD$max.iter
  } else if (meth == "CGMM") {
    isTRUE(r$fit$models[[1]]$iter.stop < COMPETITOR_DEFAULTS$CNmixt$iter.max)
  } else if (grepl("MixSIT", meth)) {
    isTRUE(r$fit$converged)
  } else NA
  ic <- if (ok) calc_information_criteria(d$X, r$fit, meth, d$K, d$m)
        else list(AIC = NA_real_, BIC = NA_real_, ICL = NA_real_)
  rows[[length(rows) + 1]] <- data.frame(
    dataset = nm, n = d$n, m = d$m, K = d$K, method = meth,
    ARI = if (ok) mclust::adjustedRandIndex(d$y, r$cluster) else NA_real_,
    NMI = if (ok) calc_nmi(d$y, r$cluster) else NA_real_,
    BIC = ic$BIC, ICL = ic$ICL, secs = secs, success = ok, converged = conv,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-12s %7.1fs  ARI %s\n", meth, secs,
              if (ok) sprintf("%.4f", tail(rows, 1)[[1]]$ARI) else "FAIL"))
  flush.console()
}

save_result(do.call(rbind, rows), out)
cat("written", out, "\n")
