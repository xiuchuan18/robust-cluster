args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
figdir <- file.path(root, "Stage2", "figures")

expected <- sort(c(
  "Figure1_ARI_Boxplots.pdf",
  "Figure2_BIC_Boxplots.pdf",
  "Section4_Wbar_Density.pdf",
  "Section5_CellWise05_ARI.pdf",
  "Section6_Scaling.pdf",
  file.path("calibrated-flexibility", c(
    "CalFlex_PhaseMap.pdf"
  ))
))
actual <- sort(list.files(figdir, recursive = TRUE, all.files = FALSE))
actual <- gsub("\\\\", "/", actual)
expected <- gsub("\\\\", "/", expected)
if (!identical(actual, expected)) {
  stop("final figure set mismatch\nexpected: ", paste(expected, collapse = ", "),
       "\nactual: ", paste(actual, collapse = ", "))
}

paths <- file.path(figdir, expected)
info <- file.info(paths)
if (any(is.na(info$size)) || any(info$size < 4000)) stop("missing or implausibly small PDF")
signatures <- vapply(paths, function(path) {
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  rawToChar(readBin(con, "raw", n = 4L))
}, character(1))
if (any(signatures != "%PDF")) stop("invalid PDF signature")

pdfinfo <- Sys.which("pdfinfo")
if (nzchar(pdfinfo)) {
  for (path in paths) {
    meta <- system2(pdfinfo, shQuote(path), stdout = TRUE, stderr = TRUE)
    pages <- grep("^Pages:", meta, value = TRUE)
    if (length(pages) != 1L || !grepl("Pages:[[:space:]]+1$", pages))
      stop("expected a one-page figure: ", path)
  }
}
cat("figure integrity passed:", length(paths), "information-designed PDFs\n")

rows <- read.csv(file.path(root,'Stage2/results/competitor-extension/simulation_plot_data.csv'))
stopifnot(nrow(rows)==7800,all(table(rows$method,rows$scenario)==100))
s <- read.csv(file.path(root,'Stage2/results/competitor-extension/section6.csv'))
stopifnot(sum(s$section=='scale')==360,sum(s$returned[s$section=='scale'])==350,
          sum(is.finite(s$angle_deg[s$section=='recovery']))==250)
