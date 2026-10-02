args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
lib <- file.path(root, ".r-lib")
if (dir.exists(lib)) .libPaths(c(lib, .Library))
declared <- utils::read.csv(file.path(root, "software-versions.csv"),
                            stringsAsFactors = FALSE)
pkgs <- setdiff(declared$package, "R")
lines <- c(capture.output(sessionInfo()), "", "Direct package versions:",
  sprintf("%-20s %s", pkgs, vapply(pkgs, function(p)
    as.character(utils::packageVersion(p)), character(1))))
out <- file.path(root, "docs", "sessionInfo.txt")
writeLines(lines, out)
cat("wrote", out, "\n")
