args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
lib <- file.path(root, ".r-lib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .Library))

required_r <- "4.4.1"
if (as.character(getRversion()) != required_r)
  stop("This release is pinned to R ", required_r, "; found ", getRversion())

cran <- "https://cloud.r-project.org"
versions <- utils::read.csv(file.path(root, "software-versions.csv"),
                            stringsAsFactors = FALSE)

installed_local <- function(package, version) {
  p <- tryCatch(find.package(package, lib.loc = lib), error = function(e) "")
  nzchar(p) && identical(as.character(utils::packageVersion(package,
    lib.loc = lib)), version)
}

# Bootstrap the installer in the local library, then install every declared
# direct dependency there. Global user libraries are deliberately ignored.
if (!installed_local("remotes", "2.5.0"))
  install.packages("remotes", lib = lib, repos = cran)
if (as.character(utils::packageVersion("remotes", lib.loc = lib)) != "2.5.0")
  remotes::install_version("remotes", "2.5.0", lib = lib, repos = cran,
                           upgrade = "never")

cran_rows <- versions[versions$source == "CRAN", , drop = FALSE]
for (i in seq_len(nrow(cran_rows))) {
  p <- cran_rows$package[i]; v <- cran_rows$version[i]
  if (!installed_local(p, v))
    remotes::install_version(p, version = v, lib = lib, repos = cran,
                             upgrade = "never")
}

if (!installed_local("EMMIXuskew", "0.11-6"))
  remotes::install_url(paste0(
    "https://cran.r-project.org/src/contrib/Archive/EMMIXuskew/",
    "EMMIXuskew_0.11-6.tar.gz"), lib = lib, upgrade = "never")

# Always build the documented compatibility patch into the local library; a
# same-version global EMMIXskew installation is not interchangeable.
install.packages(file.path(root, "vendor", "EMMIXskew"), lib = lib,
                 repos = NULL, type = "source")

old <- getwd()
setwd(file.path(root, "Stage1"))
on.exit(setwd(old), add = TRUE)
source("Build_MixSITcore_Package.R")
cat("runtime ready at", lib, "\n")
