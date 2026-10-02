# 00_RNG.R
#
# The randomisation standard for the whole study. Sourced by 00_CONFIG.R, which calls rng_fix() at load time, so every script that sources the config is
# already under the standard.


RNG_SCHEME_VERSION <- "2026-09-05.1"
EXPERIMENT_RUN_ID <- "mixsit-final-nstart5-eta5-2026-09-25"
MODEL_REVISION <- "nu-actual-baseline-2026-10-02"
CRITERIA_REVISION <- "msghd5d-softICL-posterior-matrix-2026-10-02"
PROVENANCE_PACKAGE_VERSIONS <- c(
  mclust = "6.1.1", EMMIXskew = "1.0.3", MixGHD = "2.3.7",
  ContaminatedMixt = "1.3.8", processx = "3.8.6", digest = "0.6.37")

# Pin the generator. Called by 00_CONFIG.R at load, and again by any worker
# process, so that a stream depends on its seed and on nothing else.
rng_fix <- function() {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  invisible(RNGkind())
}

mk_seed <- function(...) {
  key <- paste(as.character(unlist(list(...))), collapse = "|")
  strtoi(substr(digest::digest(key, algo = "md5", serialize = FALSE), 1, 7), 16L)
}

# Set the seed for one role and return it, so a caller can record what it used.
seed_for <- function(...) { s <- mk_seed(...); set.seed(s); s }

# Stamped onto every results object, so a file that will not reproduce can be told apart from a file that was written under different software.
rng_provenance <- function() {
  pv <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                             error = function(e) NA_character_)
  list(run_id       = EXPERIMENT_RUN_ID,
       model_revision = MODEL_REVISION,
       criteria_revision = CRITERIA_REVISION,
       scheme       = RNG_SCHEME_VERSION,
       rngkind      = RNGkind(),
       r_version    = R.version.string,
       packages     = c(mclust = pv("mclust"), EMMIXskew = pv("EMMIXskew"),
                        MixGHD = pv("MixGHD"),
                        ContaminatedMixt = pv("ContaminatedMixt"),
                        processx = pv("processx"), digest = pv("digest")),
       stamped_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
}

# Write in the destination directory and rename only after serialization has completed. A killed worker therefore cannot leave a readable partial cache.
save_result <- function(object, file, compress = TRUE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  attr(object, "provenance") <- rng_provenance()
  tmp <- tempfile(paste0(".", basename(file), "."), tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("could not atomically publish result: ", file)
  invisible(file)
}

cache_valid <- function(file) {
  if (!file.exists(file)) return(FALSE)
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(x)) return(FALSE)
  p <- attr(x, "provenance", exact = TRUE)
  affected <- grepl("MixSIT_(None|U)|^_rec_(EV|VV)_|^_lam_|^_rep_|^_real", basename(file))
  metric_affected <- grepl("MSGHD", basename(file)) && !grepl("^_scale_", basename(file))
  is.list(p) && (!affected || (identical(p$model_revision, MODEL_REVISION) &&
    identical(p$criteria_revision, CRITERIA_REVISION))) &&
    (!metric_affected || identical(p$criteria_revision, CRITERIA_REVISION)) &&
    identical(p$run_id, EXPERIMENT_RUN_ID) &&
    identical(p$scheme, RNG_SCHEME_VERSION) &&
    identical(unname(p$packages[names(PROVENANCE_PACKAGE_VERSIONS)]),
              unname(PROVENANCE_PACKAGE_VERSIONS))
}

write_csv_atomic <- function(object, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(".", basename(file), "."), tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(object, tmp, row.names = FALSE, ...)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("could not atomically publish CSV: ", file)
  invisible(file)
}
