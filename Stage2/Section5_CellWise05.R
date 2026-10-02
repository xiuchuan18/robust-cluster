# Section5_CellWise05.R
#
#   Rscript Section5_CellWise05.R [n_parallel]
#
# Runs three benchmark datasets 50 independent 5% cell-wise contamination replicates.  Each dataset-replicate is isolated in an operating-system
# process and contains all methods under the common comparison protocol. Existing per-replicate files are skipped, so the sweep resumes safely.

setwd(normalizePath(".", winslash = "/"))
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source("Section5_Datasets.R")
  source("Section5_ContaminationLib.R")
}))

args      <- commandArgs(trailingOnly = TRUE)
N_PAR     <- if (length(args) >= 1L) as.integer(args[1]) else 8L
TIMEOUT_S <- 7200
POLL_S    <- 5
if (!is.finite(N_PAR) || N_PAR < 1L) stop("n_parallel must be positive")

grid <- expand.grid(dataset = SECTION5_DATASETS,
                    rep = seq_len(SECTION5_CONTAM_NREP),
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$rep, match(grid$dataset, SECTION5_DATASETS)), ]
rp <- section5_contam_cache
done <- mapply(function(nm, r) cache_valid(rp(nm, r)), grid$dataset, grid$rep)
todo <- grid[!done, , drop = FALSE]

run_process_queue(todo, "Section5_RunContaminatedRep.R",
  args_fun = function(j) c(j$dataset, j$rep),
  done_fun = function(j) cache_valid(rp(j$dataset, j$rep)),
  label_fun = function(j) sprintf("section5_contam_%s_%02d", j$dataset, j$rep),
  n_parallel = N_PAR, timeout_s = TIMEOUT_S, poll_s = POLL_S)

files <- mapply(rp, grid$dataset, grid$rep)
missing <- files[!vapply(files, cache_valid, logical(1))]
if (length(missing))
  stop("missing Section 5 contaminated results: ",
       paste(basename(missing), collapse = ", "))
res <- do.call(rbind, lapply(files, readRDS))
save_result(res, section5_contam_result("rds"))
write_csv_atomic(res, section5_contam_result("csv"))

M <- c("GMM", "tMM", "SNMM", "STMM", "MSGHD", "CGMM",
       "MixSIT_None", "MixSIT_U", "MixSIT_Full")
cat("================= 5% cell-wise contamination =================\n")
cat(sprintf("entries: mean ARI (SD), with successful fits / %d\n\n",
            SECTION5_CONTAM_NREP))
for (nm in SECTION5_DATASETS) {
  cat(nm, "\n")
  for (meth in M) {
    z <- res[res$dataset == nm & res$method == meth, ]
    cat(sprintf("  %-12s %7.3f (%5.3f)  %2d/%d\n",
                METHOD_CODE_TO_LABEL[[meth]], mean(z$ARI, na.rm = TRUE),
                stats::sd(z$ARI, na.rm = TRUE), sum(is.finite(z$ARI)), nrow(z)))
  }
  u <- unique(res[res$dataset == nm,
                  c("rep", "fraction_contaminated", "fraction_rows_affected")])
  cat(sprintf("  realised contamination: %.3f of cells; %.3f of rows affected\n\n",
              mean(u$fraction_contaminated), mean(u$fraction_rows_affected)))
}
