# Shared definitions for the confirmatory calibrated-flexibility experiments.
#
# These experiments use the production estimator and competitor wrappers from
# Section4_Lib.R, but have a distinct run id and seed namespace so exploratory
# caches can never be promoted accidentally.

CALFLEX_RUN_ID <- "mixsit-calflex-formal-2026-09-26"
CALFLEX_SEED_PREFIX <- "CalFlexFormal20260926"

calflex_provenance <- function() {
  p <- rng_provenance()
  p$run_id <- CALFLEX_RUN_ID
  p$parent_run_id <- EXPERIMENT_RUN_ID
  p
}

calflex_save <- function(object, file, compress = TRUE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  attr(object, "provenance") <- calflex_provenance()
  tmp <- tempfile(paste0(".", basename(file), "."), tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("could not publish result: ", file)
  invisible(file)
}

calflex_cache_valid <- function(file) {
  if (!file.exists(file)) return(FALSE)
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(x)) return(FALSE)
  p <- attr(x, "provenance", exact = TRUE)
  affected <- grepl("MixSIT_(None|U)|banknote|Banknote", basename(file))
  is.list(p) && (!affected || (identical(p$model_revision, MODEL_REVISION) &&
    identical(p$criteria_revision, CRITERIA_REVISION))) &&
    (!grepl("MSGHD", basename(file)) || identical(p$criteria_revision, CRITERIA_REVISION)) &&
    identical(p$run_id, CALFLEX_RUN_ID) &&
    identical(p$scheme, RNG_SCHEME_VERSION)
}

calflex_write_csv <- function(object, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(".", basename(file), "."), tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(object, tmp, row.names = FALSE, ...)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("could not publish CSV: ", file)
  invisible(file)
}

calflex_seed <- function(...) mk_seed(CALFLEX_SEED_PREFIX, ...)
calflex_fit_seed <- function(tag, method, rep) mk_seed(tag, "fit", method, rep)

calflex_methods <- c("GMM", "tMM", "SNMM", "STMM", "MSGHD", "CGMM",
                     "MixSIT_None", "MixSIT_U", "MixSIT_Full")
calflex_phase_methods <- c("GMM", "MSGHD", "MixSIT_None", "MixSIT_U",
                           "MixSIT_Full")

calflex_method_label <- METHOD_CODE_TO_LABEL

calflex_load_data <- function(name) {
  if (name %in% SECTION5_DATASETS) return(load_dataset(name))
  if (name != "Banknote") stop("unknown calibrated-flexibility dataset: ", name)
  e <- new.env(parent = emptyenv())
  utils::data("banknote", package = "mclust", envir = e)
  b <- e$banknote
  X <- scale(as.matrix(b[, setdiff(names(b), "Status"), drop = FALSE]))
  y <- as.integer(b$Status)
  list(name = name, X = X, y = y, K = nlevels(b$Status), n = nrow(X),
       m = ncol(X), n_duplicate = sum(duplicated(X)), dedup = FALSE,
       why = "six measurements from the Swiss banknote data")
}

calflex_converged <- function(res, method) {
  if (is.null(res) || is.null(res$fit)) return(FALSE)
  f <- res$fit
  if (method %in% c("tMM", "SNMM", "STMM")) isTRUE(f$error == 0)
  else if (method == "MSGHD")
    length(f@loglik) < COMPETITOR_DEFAULTS$MSGHD$max.iter
  else if (method == "CGMM")
    isTRUE(f$models[[1]]$iter.stop < COMPETITOR_DEFAULTS$CNmixt$iter.max)
  else if (grepl("MixSIT", method)) isTRUE(f$converged)
  else NA
}

# MSGHD's fitted component contains mu, alpha, phi, omega and lambda (five d-vectors), plus an orthogonal Gamma matrix.  MixGHD 2.3.7 counts only 4d
# and omits phi.  An earlier exploratory audit used 6d; that over-corrected by one d-vector.  The formal primary criterion uses the actual 5d count and
# retains both alternatives as sensitivity columns.
calflex_msghd_npar <- function(K, d, convention = c("actual", "package", "over")) {
  convention <- match.arg(convention)
  nd <- switch(convention, actual = 5L, package = 4L, over = 6L)
  (K - 1) + K * (nd * d + d * (d - 1) / 2)
}

calflex_npar <- function(method, K, d) {
  count_parameters(method, K, d)
}

calflex_entropy <- function(z) posterior_entropy(z)

calflex_posterior <- function(fit, method) {
  model_posterior(fit, method)
}
