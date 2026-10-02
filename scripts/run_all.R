args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
lib <- file.path(root, ".r-lib")
if (dir.exists(lib)) .libPaths(c(lib, .Library))
Sys.setenv(R_LIBS_USER = lib, OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
stage2 <- file.path(root, "Stage2")
setwd(stage2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
logfile <- file.path("results", "run_all.log")
logcon <- file(logfile, open = "at")
sink(logcon, type = "output")
sink(logcon, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(logcon)}, add = TRUE)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")
stamp <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           sprintf(...), "\n")
run <- function(script, arguments = character()) {
  stamp("START %s %s", script, paste(arguments, collapse = " "))
  t0 <- Sys.time()
  status <- system2(rscript, c(script, arguments))
  stamp("END   %s (%.1f min; status %d)", script,
        as.numeric(difftime(Sys.time(), t0, units = "mins")), status)
  if (status != 0L) stop(script, " failed with status ", status)
}

stamp("FORMAL EXPERIMENT WORKFLOW START")
run(file.path("..", "tests", "preflight.R"))
run(file.path("..", "tests", "method_revision.R"))
run(file.path("..", "tests", "information_criteria.R"))
run(file.path("..", "tests", "wbar_axis_matching.R"))
run(file.path("..", "tests", "seed_audit.R"))
run(file.path("..", "tests", "calflex_seed_audit.R"))
# Parallelism is staged by memory footprint.  These experiments use semantic
# per-job RNG streams, so scheduling cannot change their numerical draws.  Only
# the scalability study below is a timing experiment and remains fixed at eight
# concurrent fits, exactly as stated in the manuscript.
run("Section4_Run.R", "8")
run("Section4_BayesCeiling.R")
run("Section5_RealData.R", "3")
run("Section5_CellWise05.R", "12")
run("Section6_Recovery.R", "12")
run("Section6_LambdaSensitivity.R", "12")
run("Section6_KSelect.R", "10")

# Timing study runs alone under a fixed eight-process load.
run("Section6_Scale.R", "8")

# Confirmatory calibrated-flexibility experiments are part of the same formal
# study.  They run here, after the Section 5 data products they reuse and before
# the combined figures are assembled, so a clean checkout has no circular
# dependency between the two historical runners.
run("CalFlex_BayesProfile.R")
run("CalFlex_Banknote.R", c("driver", "14"))
run("CalFlex_PhaseMap.R", c("driver", "14"))

for (script in c("Section4_Figures.R", "Section4_WbarAnalysis.R",
                 "CalFlex_BanknoteAnalysis.R", "CalFlex_PhaseMapAnalysis.R",
                 "Section5_CellWise05_Figure.R", "Final_Summaries.R",
                 "Section6_ScalingFigure.R")) run(script)

run(file.path("..", "tests", "result_integrity.R"))
run(file.path("..", "tests", "calflex_integrity.R"))
run(file.path("..", "tests", "figure_integrity.R"))
run(file.path("..", "scripts", "write_session_info.R"))
run(file.path("..", "scripts", "write_results_manifest.R"))
stamp("FORMAL EXPERIMENT WORKFLOW COMPLETE")
