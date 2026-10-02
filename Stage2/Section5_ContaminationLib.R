# Section5_ContaminationLib.R
#
# Cell-wise contamination used for the semi-real extension of Section 5. The clean benchmark variables are already standardized by load_dataset().
# Each standardized cell is independently contaminated with probability alpha by adding N(0, eta^2) noise.  The data are not re-standardized afterwards:
# doing so would let the contaminated sample alter the scale definition and would make eta incomparable with the controlled Section 4 experiment.

SECTION5_CONTAM_ALPHA <- 0.05
SECTION5_CONTAM_ETA   <- 5
SECTION5_CONTAM_NREP  <- 50L
SECTION5_CONTAM_ID    <- "cellwise05_eta5"
SECTION5_CONTAM_STEM  <- "Section5_CellWise05"

section5_contam_cache <- function(dataset, rep) {
  file.path(RESULTS_DIR,
            sprintf("_realcontam_a05_e5_%s_%02d.rds", dataset, rep))
}

section5_contam_result <- function(ext = "rds") {
  file.path(RESULTS_DIR, paste0(SECTION5_CONTAM_STEM, ".", ext))
}

contaminate_section5 <- function(X, dataset, rep,
                                 alpha = SECTION5_CONTAM_ALPHA,
                                 eta = SECTION5_CONTAM_ETA) {
  stopifnot(is.matrix(X), length(dataset) == 1L, length(rep) == 1L,
            is.finite(alpha), alpha >= 0, alpha <= 1,
            is.finite(eta), eta >= 0)
  seed <- seed_for("Section5", dataset, "cellwise", alpha, eta, "contam", rep)
  mask <- matrix(stats::runif(length(X)) < alpha, nrow(X), ncol(X))
  shock <- matrix(stats::rnorm(length(X), sd = eta), nrow(X), ncol(X))
  Xc <- X + mask * shock
  list(
    X = Xc,
    mask = mask,
    seed = seed,
    n_contaminated = sum(mask),
    fraction_contaminated = mean(mask),
    n_rows_affected = sum(rowSums(mask) > 0),
    fraction_rows_affected = mean(rowSums(mask) > 0)
  )
}
