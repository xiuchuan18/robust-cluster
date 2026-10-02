# Section5_RealData.R
#
#   Rscript Section5_RealData.R [n_parallel]
#
# Drives Section5_RunDataset.R over the three formal, unmodified benchmark datasets (Biopsy, full 11-variable AIS, Wine), one operating-system process
# each, and assembles results/Section5_RealData.rds. Re-running skips datasets that already have a result file, so an interrupted
# sweep resumes.

setwd(normalizePath(".", winslash = "/"))
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source("Section5_Datasets.R")
}))

args      <- commandArgs(trailingOnly = TRUE)
N_PAR     <- if (length(args) >= 1) as.integer(args[1]) else 6
TIMEOUT_S <- 3600
POLL_S    <- 5

NAMES <- SECTION5_DATASETS
rp <- function(nm) file.path(RESULTS_DIR, paste0("_real_", nm, ".rds"))
todo <- NAMES[!vapply(NAMES, function(nm) cache_valid(rp(nm)), logical(1))]

jobs <- data.frame(dataset = todo, stringsAsFactors = FALSE)
run_process_queue(jobs, "Section5_RunDataset.R",
  args_fun = function(j) j$dataset,
  done_fun = function(j) cache_valid(rp(j$dataset)),
  label_fun = function(j) paste0("section5_clean_", j$dataset),
  n_parallel = N_PAR, timeout_s = TIMEOUT_S, poll_s = POLL_S)

missing <- NAMES[!vapply(NAMES, function(nm) cache_valid(rp(nm)), logical(1))]
if (length(missing)) stop("missing Section 5 result files: ", paste(missing, collapse = ", "))
res <- do.call(rbind, lapply(NAMES, function(nm) readRDS(rp(nm))))
save_result(res, file.path(RESULTS_DIR, "Section5_RealData.rds"))
write_csv_atomic(res, file.path(RESULTS_DIR, "Section5_RealData.csv"))

# report
M <- c("GMM", "tMM", "SNMM", "STMM", "MSGHD", "CGMM",
       "MixSIT_None", "MixSIT_U", "MixSIT_Full")
L <- vapply(M, function(m) METHOD_CODE_TO_LABEL[[m]], "")
cat("=========================== ARI ===========================\n")
cat(sprintf("%-9s %4s %3s %2s  %s\n", "dataset", "n", "m", "K",
            paste(sprintf("%9s", L), collapse = "")))
for (nm in unique(res$dataset)) {
  d <- res[res$dataset == nm, ]
  a <- vapply(M, function(m) d$ARI[d$method == m], 0)
  cat(sprintf("%-9s %4d %3d %2d  %s\n", nm, d$n[1], d$m[1], d$K[1],
              paste(vapply(a, function(v)
                if (is.na(v)) sprintf("%9s", "--") else sprintf("%9.4f", v), ""),
                collapse = "")))
}

cat("\n============== best MixSIT vs best competitor ==============\n")
cat(sprintf("%-9s %14s %-12s %14s %-10s %9s\n",
            "dataset", "best MixSIT", "", "best compet.", "", "margin"))
for (nm in unique(res$dataset)) {
  d <- res[res$dataset == nm, ]
  a <- vapply(M, function(m) d$ARI[d$method == m], 0); names(a) <- M
  a[is.na(a)] <- -Inf
  mx <- a[grepl("MixSIT", names(a))]; cp <- a[!grepl("MixSIT", names(a))]
  cat(sprintf("%-9s %14.4f %-12s %14.4f %-10s %+9.4f\n", nm,
              max(mx), METHOD_CODE_TO_LABEL[[names(which.max(mx))]],
              max(cp), METHOD_CODE_TO_LABEL[[names(which.max(cp))]],
              max(mx) - max(cp)))
}

cat("\n=================== seconds per fit ===================\n")
cat(sprintf("%-9s %s\n", "dataset", paste(sprintf("%9s", L), collapse = "")))
for (nm in unique(res$dataset)) {
  d <- res[res$dataset == nm, ]
  cat(sprintf("%-9s %s\n", nm,
      paste(vapply(M, function(m) sprintf("%9.1f", d$secs[d$method == m]), ""),
            collapse = "")))
}
