# 00_CONFIG.R
#
# Global configuration shared by all Stage 2 experiment scripts.
# Centralizes seeds, lambda values, max_iter, paths, and parallel setup.
# Source this file at the top of every experiment script.

# Prefer the repository-local runtime library when present. It is excluded from
# version control; scripts/setup_runtime.R recreates it from pinned versions.
LOCAL_R_LIB <- normalizePath(file.path("..", ".r-lib"), winslash = "/",
                             mustWork = FALSE)
if (dir.exists(LOCAL_R_LIB)) .libPaths(c(LOCAL_R_LIB, .Library))

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(digest)
})

# Randomisation standard
# Defines mk_seed(), seed_for(), rng_fix(), save_result(). rng_fix() is called immediately so every script that sources this config is already under the
# pinned generator.
source("00_RNG.R")   # every script sources this config from Stage2/
rng_fix()
source("00_ProcessQueue.R")

# Paths
STAGE1_DIR <- "../Stage1"
RESULTS_DIR <- "results"
FIGURES_DIR <- "figures"

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(FIGURES_DIR)) dir.create(FIGURES_DIR, recursive = TRUE)

# MixSIT model settings
MIXSIT <- list(
  max_iter         = 500,
  tol              = 1e-4,
  n_short          = 10,
  short_iter       = 25,
  init_methods     = c("kmeans", "hclust", "random"),
  lambda_U         = 0,
  lambda_O         = 0,
  lambda_magnitude = 0,
  lambda_hetero    = 5,
  variants         = c("none", "U", "full"),   # internal codes
  variant_labels   = c("none" = "MixSIT-VV",
                       "U"    = "MixSIT-EV",
                       "full" = "MixSIT-EE")
)

# Cross-section competitor protocol
# One rule, applied identically in Sections 4, 5 and 6. It is stated in three parts because a single "use the package default" would be ambiguous where a
# package's default is itself a multi-start.

COMPETITOR_MAX_ITER <- 500      # MixSIT's cap; competitors use their own defaults
COMPETITOR_TOL      <- 1e-4     # ceiling, not target: see TOLERANCE above
KMEANS_NSTART       <- 5L

# Verified package defaults, kept here so a version bump that changes one is
# caught by COMPETITOR_DEFAULTS_CHECK() rather than silently altering the study.
COMPETITOR_DEFAULTS <- list(
  MSGHD  = list(method = "km", max.iter = 100, eps = 0.01, modelSel = "AIC",
                scale = TRUE),
  EmSkew = list(nkmeans = 0, nrandom = 10, nhclust = FALSE, initloop = 20,
                itmax = 1000, epsilon = 1e-6),
  CNmixt = list(initialization = "mixt", iter.max = 1000)
)

# Re-read the installed packages and warn if any default has moved.
COMPETITOR_DEFAULTS_CHECK <- function() {
  chk <- function(pkg, fname, want) {
    if (!requireNamespace(pkg, quietly = TRUE)) return("package not installed")
    a <- tryCatch(formals(getExportedValue(pkg, fname)), error = function(e) NULL)
    if (is.null(a)) return("no formals")
    bad <- names(want)[vapply(names(want),
             function(k) !isTRUE(all.equal(a[[k]], want[[k]])), TRUE)]
    if (length(bad)) paste0("CHANGED (", paste(bad, collapse = ", "), ")") else "ok"
  }
  c(MSGHD  = chk("MixGHD", "MSGHD", COMPETITOR_DEFAULTS$MSGHD),
    EmSkew = chk("EMMIXskew", "EmSkew", COMPETITOR_DEFAULTS$EmSkew),
    CNmixt = chk("ContaminatedMixt", "CNmixt", COMPETITOR_DEFAULTS$CNmixt))
}

# Random seeds (preserved from original scripts)
# Section 4 simulation uses seed = 2026 (from Simulation_week12.Rmd)
# Section 6 experiments use seed = 20260208 (from Combined_5_Experiments.R)
# Each experiment derives per-replicate seeds by adding fixed offsets.
SEED_SIM    <- 2026
SEED_MECH   <- 20260208

# Parallel settings
N_CORES <- max(1, parallel::detectCores() - 2)

# Method labels (manuscript-consistent)
# Maps internal code names -> labels used in manuscript figures and tables.
METHOD_CODE_TO_LABEL <- c(
  "GMM"          = "GMM",
  "tMM"          = "tMM",
  "SNMM"         = "MSN",
  "STMM"         = "MST",
  "MSGHD"        = "MSGHD",
  "CGMM"         = "CGMM",
  "MixSIT_None"  = "MixSIT-VV",
  "MixSIT_U"     = "MixSIT-EV",
  "MixSIT_Full"  = "MixSIT-EE"
)

METHOD_LABEL_ORDER <- c(
  "GMM", "tMM", "MSN", "MST", "MSGHD", "CGMM",
  "MixSIT-VV", "MixSIT-EV", "MixSIT-EE"
)

# Visualization colors
METHOD_COLORS <- c(
  "GMM"        = "#4E79A7",
  "tMM"        = "#76B7B2",
  "MSN"        = "#F28E2B",
  "MST"        = "#E15759",
  "MSGHD"      = "#59A14F",
  "CGMM"       = "#9C755F",
  "MixSIT-VV"  = "#1B9E77",
  "MixSIT-EV"  = "#66C2A5",
  "MixSIT-EE"  = "#A6D854"
)

# Helper: source Stage 1 core code
# Each experiment script must call this at startup so MixSIT functions are available in the master R process AND in parallel workers.
source_stage1 <- function() {
  prev_wd <- getwd()
  setwd(STAGE1_DIR)
  source("MixSIT_rcpp.R")
  source("MixSIT_flip.R")   # the axis-flip stage is part of the estimator
  setwd(prev_wd)
}

# Helper: convenience ARI
# Unlabelled observations are dropped from both vectors rather than passed through. Every model-based arm labels the whole sample, so for those this is a
# no-op -- their cluster vector is either complete or entirely NA. It matters only for the DDC_MNM two-step arm, where DDC can strip a row of all but one
# entry and the mixture then has nothing to classify it from. Scoring such a row as a miss would penalise the baseline for a decision the pipeline is entitled
# to make; leaving NA in place would make adjustedRandIndex treat it as its own cluster, which is worse. The count of dropped rows is recorded alongside the
# ARI so the size of the subset is always visible.
.drop_unlabelled <- function(true_labels, pred_labels) {
  keep <- !is.na(pred_labels)
  list(t = true_labels[keep], p = pred_labels[keep], n = sum(keep))
}

calc_ari <- function(true_labels, pred_labels) {
  if (is.null(pred_labels) || length(pred_labels) != length(true_labels)) return(NA_real_)
  if (all(is.na(pred_labels))) return(NA_real_)
  z <- .drop_unlabelled(true_labels, pred_labels)
  if (z$n < 10) return(NA_real_)
  mclust::adjustedRandIndex(z$t, z$p)
}

calc_nmi <- function(true_labels, pred_labels) {
  if (is.null(pred_labels) || length(pred_labels) != length(true_labels)) return(NA_real_)
  if (all(is.na(pred_labels))) return(NA_real_)
  z <- .drop_unlabelled(true_labels, pred_labels)
  if (z$n < 10) return(NA_real_)
  tryCatch(aricode::NMI(z$t, z$p), error = function(e) NA_real_)
}

cat("CONFIG loaded: max_iter=", MIXSIT$max_iter,
    "  n_cores=", N_CORES,
    "  seeds: sim=", SEED_SIM, " mech=", SEED_MECH, "\n", sep = "")
