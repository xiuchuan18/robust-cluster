args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
lib <- file.path(root, ".r-lib")
if (dir.exists(lib)) .libPaths(c(lib, .Library))
setwd(file.path(root, "Stage2"))

if (as.character(getRversion()) != "4.4.1")
  stop("R 4.4.1 is required; found ", getRversion())

declared <- utils::read.csv(file.path(root, "software-versions.csv"),
                            stringsAsFactors = FALSE)
declared <- declared[!declared$package %in% "R", , drop = FALSE]
for (i in seq_len(nrow(declared))) {
  p <- declared$package[i]; want <- declared$version[i]
  if (!requireNamespace(p, quietly = TRUE)) stop("missing package: ", p)
  got <- as.character(utils::packageVersion(p))
  if (!identical(package_version(got), package_version(want)))
    stop(p, " version mismatch: expected ", want, ", found ", got)
}

suppressWarnings(suppressMessages({
  source("00_CONFIG.R")
  source_stage1()
  source("Section4_Lib.R")
}))

stopifnot(KMEANS_NSTART == 5L, R_RESTART == 1L)
defaults <- COMPETITOR_DEFAULTS_CHECK()
if (!all(defaults == "ok"))
  stop("competitor defaults changed: ", paste(names(defaults), defaults,
                                               collapse = "; "))

# Verify that the vendored source and its checksum manifest agree.
vendor <- file.path(root, "vendor", "EMMIXskew")
md5_lines <- readLines(file.path(vendor, "MD5"))
md5_expected <- sub(" .*", "", md5_lines)
md5_files <- sub("^[^ ]+ \\*", "", md5_lines)
md5_actual <- unname(tools::md5sum(file.path(vendor, md5_files)))
if (!identical(md5_actual, md5_expected)) stop("vendored EMMIXskew MD5 mismatch")

# nstart=1 must exactly reproduce the native MSGHD method='km' initializer.
Y <- scale(as.matrix(iris[1:100, 1:4]))
set.seed(77)
native <- .quiet(MixGHD::MSGHD(data = Y, G = 2, method = "km",
  max.iter = 2, eps = 1e-4, scale = FALSE, modelSel = "ICL"))
set.seed(77)
gpar <- .msghd_kmeans_gpar(Y, 2, nstart = 1L)
wrapped <- .quiet(MixGHD::MSGHD(data = Y, gpar0 = gpar, G = 2,
  method = "km", max.iter = 2, eps = 1e-4, scale = FALSE, modelSel = "ICL"))
stopifnot(identical(native@map, wrapped@map))
stopifnot(max(abs(native@loglik - wrapped@loglik)) == 0)

for (distr in c("mvt", "msn", "mst")) {
  set.seed(99)
  ini <- .emmix_kmeans_init(Y, 2, distr, nstart = 5L)
  stopifnot(!is.null(ini))
  fit <- .quiet(EMMIXskew::EmSkew(Y, g = 2, distr = distr, ncov = 3,
    init = ini, itmax = 3, epsilon = 1e-6, debug = FALSE))
  stopifnot(!is.null(fit))
}

source("Section5_Datasets.R")
expected_dims <- c(Biopsy = 9L, AIS = 11L, Wine = 13L)
for (name in names(expected_dims)) {
  dat <- load_dataset(name)
  stopifnot(ncol(dat$X) == expected_dims[[name]], nrow(dat$X) == length(dat$y),
            all(is.finite(dat$X)))
}

cat("preflight passed\n")
