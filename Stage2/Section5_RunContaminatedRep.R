# Section5_RunContaminatedRep.R
#
#   Rscript Section5_RunContaminatedRep.R <dataset-name> <replicate>
#
# Fits all methods to one independently contaminated copy of one benchmark dataset.  Initialization, iteration caps, tolerances and failure handling are
# exactly those used by Sections 4 and 5.  One file is written atomically per dataset-replicate so the 90-job sweep is resumable.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: Rscript Section5_RunContaminatedRep.R <dataset> <replicate>")
nm  <- args[1]
rep <- as.integer(args[2])
if (!is.finite(rep) || rep < 1L) stop("replicate must be a positive integer")

setwd(normalizePath(".", winslash = "/"))

suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R")
  source("Section5_Datasets.R"); source("Section5_ContaminationLib.R")
}))

out <- section5_contam_cache(nm, rep)

d <- load_dataset(nm)
cc <- contaminate_section5(d$X, nm, rep)
cat(sprintf("%s rep %02d: n=%d m=%d K=%d, %d/%d cells (%.3f), %d/%d rows (%.3f)\n",
            nm, rep, d$n, d$m, d$K, cc$n_contaminated, length(cc$mask),
            cc$fraction_contaminated, cc$n_rows_affected, d$n,
            cc$fraction_rows_affected))
flush.console()

rows <- list()
rng_fix()
for (meth in METHOD_NAMES) {
  t0 <- Sys.time()
  fit <- tryCatch(
    fit_multistart(meth, cc$X, d$K, R_RESTART, COMPETITOR_MAX_ITER,
                   tag = paste("Section5CellWise05", nm, sep = "/"), rep = rep),
    error = function(e) structure(list(message = conditionMessage(e)),
                                  class = "section5_fit_error"))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  errored <- inherits(fit, "section5_fit_error")
  ok <- !errored && !is.null(fit) && !is.null(fit$cluster) &&
        !all(is.na(fit$cluster))
  conv <- if (!ok) FALSE else if (meth %in% c("tMM", "SNMM", "STMM")) {
    isTRUE(fit$fit$error == 0)
  } else if (meth == "MSGHD") {
    length(fit$fit@loglik) < COMPETITOR_DEFAULTS$MSGHD$max.iter
  } else if (meth == "CGMM") {
    isTRUE(fit$fit$models[[1]]$iter.stop < COMPETITOR_DEFAULTS$CNmixt$iter.max)
  } else if (grepl("MixSIT", meth)) {
    isTRUE(fit$fit$converged)
  } else NA
  ic <- if (ok) calc_information_criteria(cc$X, fit$fit, meth, d$K, d$m)
        else list(AIC = NA_real_, BIC = NA_real_, ICL = NA_real_)
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = nm, setting = SECTION5_CONTAM_ID, rep = rep,
    alpha = SECTION5_CONTAM_ALPHA, eta = SECTION5_CONTAM_ETA,
    contamination_seed = cc$seed,
    n = d$n, m = d$m, K = d$K,
    n_contaminated = cc$n_contaminated,
    fraction_contaminated = cc$fraction_contaminated,
    n_rows_affected = cc$n_rows_affected,
    fraction_rows_affected = cc$fraction_rows_affected,
    method = meth,
    ARI = if (ok) mclust::adjustedRandIndex(d$y, fit$cluster) else NA_real_,
    NMI = if (ok) calc_nmi(d$y, fit$cluster) else NA_real_,
    AIC = ic$AIC, BIC = ic$BIC, ICL = ic$ICL,
    secs = secs, success = ok, converged = conv,
    error = if (errored) fit$message else "",
    stringsAsFactors = FALSE)
  cat(sprintf("  %-12s %7.1fs  ARI %s%s\n", meth, secs,
              if (ok) sprintf("%.4f", tail(rows, 1)[[1]]$ARI) else "FAIL",
              if (ok && !isTRUE(conv) && !is.na(conv)) " [not converged]" else ""))
  flush.console()
}

save_result(do.call(rbind, rows), out)
cat("written", out, "\n")
