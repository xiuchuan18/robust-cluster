# Section5_Datasets.R
# The three clean benchmark datasets used by both Section 5 arms.  The clean arm uses these matrices unchanged after standardization; the semi-real arm
# applies contamination in Section5_ContaminationLib.R.

suppressPackageStartupMessages({
  library(sn)
  library(gclus)
  library(MASS)
})

.get_data <- function(pkg, name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = pkg, envir = env)
  get(name, envir = env)
}

DATASETS <- list(
  Biopsy = function() {
    d <- .get_data("MASS", "biopsy")
    d <- d[stats::complete.cases(d), ]
    list(X = scale(as.matrix(d[, 2:10])), y = as.integer(d$class), K = 2L,
         why = "nine cytological scores on the original biopsy data")
  },
  AIS = function() {
    d <- .get_data("sn", "ais")
    vars <- c("RCC", "WCC", "Hc", "Hg", "Fe", "BMI", "SSF", "Bfat",
              "LBM", "Ht", "Wt")
    list(X = scale(as.matrix(d[, vars])), y = as.numeric(d$sex), K = 2L,
         why = "all eleven physiological measurements in the original AIS data")
  },
  Wine = function() {
    d <- .get_data("gclus", "wine")
    list(X = scale(as.matrix(d[, -1])), y = as.numeric(d$Class), K = 3L,
         why = "all thirteen chemical variables in the original wine data")
  }
)

SECTION5_DATASETS <- c("Biopsy", "AIS", "Wine")

load_dataset <- function(name) {
  loader <- DATASETS[[name]]
  if (is.null(loader)) stop("unknown Section 5 dataset: ", name)
  d <- loader()
  keep <- stats::complete.cases(d$X) & !is.na(d$y)
  d$X <- d$X[keep, , drop = FALSE]
  d$y <- d$y[keep]

  # Repeated profiles in Biopsy are distinct subjects measured on a coarse integer scale; no observations are deduplicated or jittered.
  d$n_duplicate <- sum(duplicated(d$X))
  d$dedup <- FALSE
  d$name <- name
  d$n <- nrow(d$X)
  d$m <- ncol(d$X)
  d
}
