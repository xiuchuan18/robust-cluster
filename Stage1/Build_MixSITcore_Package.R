# Build_MixSITcore_Package.R
#
# Compiles a MixSIT C++ backend into an installed R package.
#
# MixSIT_rcpp.R requires the backend to be a real package rather than a sourceCpp() call, because sourceCpp() is not reliable inside parallel workers
# on Windows (each worker races to load the same temporary DLL). This script was referenced by MixSIT_rcpp.R but missing from the repository; it is restored
# here so the pipeline is reproducible from source.
#
# Usage (from the Stage1 directory):
#
#   Rscript Build_MixSITcore_Package.R
#   Rscript Build_MixSITcore_Package.R SIT_core.cpp MixSITcore
#   Rscript Build_MixSITcore_Package.R SIT_core_optimized.cpp MixSITfast
#
# The second form is what the verification script uses: it installs the optimized backend under a separate name so the two can be compared
# side by side before anything is overwritten.

args <- commandArgs(trailingOnly = TRUE)
SRC  <- if (length(args) >= 1) args[1] else "SIT_core.cpp"
PKG  <- if (length(args) >= 2) args[2] else "MixSITcore"

stopifnot(file.exists(SRC))
for (p in c("Rcpp", "RcppArmadillo")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is required to build the backend.")
}

build_dir <- file.path(tempdir(), paste0("build_", PKG, "_", as.integer(Sys.time())))
pkg_dir   <- file.path(build_dir, PKG)
dir.create(file.path(pkg_dir, "src"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(pkg_dir, "R"),   recursive = TRUE, showWarnings = FALSE)

file.copy(SRC, file.path(pkg_dir, "src", "SIT_core.cpp"), overwrite = TRUE)

# RcppArmadillo needs the R-provided LAPACK/BLAS to be linked explicitly, otherwise the solve()/inv() paths fail at link time with undefined
# references to dgetrf_, dtrtrs_, dposv_, ... -O3 over R.s default -O2. No fast-math: the floating-point semantics must not
# change, only the amount of scheduling and unrolling the compiler does.
mk <- c("PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) -O3",
        "PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)")
writeLines(mk, file.path(pkg_dir, "src", "Makevars"))
writeLines(mk, file.path(pkg_dir, "src", "Makevars.win"))

writeLines(c(
  paste0("Package: ", PKG),
  "Type: Package",
  paste0("Title: MixSIT C++ Backend (", PKG, ")"),
  "Version: 1.1.0",
  "Date: 2026-08-28",
  "Author: MixSIT authors",
  "Maintainer: MixSIT authors <noreply@example.com>",
  "Description: Compiled kernels for the MixSIT mixture model: the skew-",
  "    IStudent-t density, the E-step expectations, and the CM-step objective",
  "    functions for the three structural variants.",
  "License: GPL-3",
  "Imports: Rcpp (>= 1.0.0)",
  "LinkingTo: Rcpp, RcppArmadillo",
  "Encoding: UTF-8",
  "RoxygenNote: 7.0.0"
), file.path(pkg_dir, "DESCRIPTION"))

writeLines(c(
  paste0("useDynLib(", PKG, ", .registration=TRUE)"),
  "importFrom(Rcpp, evalCpp)",
  "exportPattern(\"^[[:alpha:]]\")"
), file.path(pkg_dir, "NAMESPACE"))

Rcpp::compileAttributes(pkg_dir)

install.packages(pkg_dir, repos = NULL, type = "source",
                 INSTALL_opts = "--no-multiarch")

if (!requireNamespace(PKG, quietly = TRUE))
  stop("Build finished but '", PKG, "' is not loadable.")

fns <- c("dskwet_cpp", "dskwet_vec_cpp", "rskwet_cpp", "rskwet_matrix_cpp",
         "vector_to_orthogonal_cpp", "E_Step_cpp", "obj_space_hetero_cpp",
         "obj_space_full_cpp", "obj_space_u_only_cpp", "obj_gamma_cpp",
         "obj_grad_hetero_cpp", "obj_grad_full_cpp", "obj_grad_u_only_cpp",
         "obj_grad_gamma_cpp", "cm2_coord_cpp")
missing <- fns[!vapply(fns, function(f) exists(f, envir = asNamespace(PKG)),
                       logical(1))]
if (length(missing))
  stop("Installed but these entry points are missing: ", paste(missing, collapse = ", "))

message("OK: '", PKG, "' installed with all ", length(fns), " entry points.")
