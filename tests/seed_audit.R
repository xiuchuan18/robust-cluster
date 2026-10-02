args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
lib <- file.path(root, ".r-lib")
if (dir.exists(lib)) .libPaths(c(lib, .Library))
setwd(file.path(root, "Stage2"))
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source("Section4_Lib.R")
  source("Section6_KSelectLib.R")
}))

keys <- character()
add <- function(...) {
  keys <<- c(keys, paste(as.character(unlist(list(...))), collapse = "|"))
}
methods <- METHOD_NAMES

for (cell in names(CONFIG_SIM$scenarios)) for (r in 1:100) {
  add(SEED_SIM, cell, "data", r)
  ds <- mk_seed(SEED_SIM, cell, "data", r)
  if (cell == "S7_CellWise10") { add(ds, "contam"); add(ds, "perm") }
  for (m in methods) add(paste(SEED_SIM, cell, sep = "/"), "fit", m, r)
}

for (nm in c("Biopsy", "AIS", "Wine")) {
  for (m in methods) add(paste("Section5", nm, sep = "/"), "fit", m, 1)
  for (r in 1:50) {
    add("Section5", nm, "cellwise", 0.05, 5, "contam", r)
    for (m in methods)
      add(paste("Section5CellWise05", nm, sep = "/"), "fit", m, r)
  }
}

for (r in 1:50) {
  add("Recovery", "D3", "truth", r)
  for (N in c(200, 400, 800, 1600, 3200)) {
    add("Recovery", "D3", "data", N, r); add("Recovery", "D3", "perm", N, r)
    for (m in c("MixSIT_U", "MixSIT_None", "MixSIT_Full"))
      add("Recovery/D3", "fit", m, paste(N, r, sep = "_"))
  }
}
for (r in 1:10) for (N in c(200, 400, 800, 1600))
  for (lam in c(0, 0.5, 1, 2, 5, 10, 25)) add("Lambda", "fit", lam, N, r)

for (cell in KCELLS) for (r in 1:30) {
  add("KSelect", cell, "data", r); add("KSelect", cell, "perm", r)
  for (m in c("GMM", "MSGHD", "MixSIT_None", "MixSIT_U", "MixSIT_Full"))
    for (K in 1:6) add(paste("KSelect", cell, sep = "/"), "fit", m,
                        paste(r, "K", K, sep = "_"))
}
cfg <- unique(rbind(
  data.frame(N = c(200, 400, 800, 1600), m = 10, K = 3),
  data.frame(N = 800, m = c(3, 5, 10, 15, 20), K = 3),
  data.frame(N = 800, m = 10, K = c(2, 3, 4, 5, 6))))
for (i in seq_len(nrow(cfg))) for (r in 1:30) {
  add("Scale", "data", cfg$N[i], cfg$m[i], cfg$K[i], r)
  tag <- paste("Scale", cfg$N[i], cfg$m[i], cfg$K[i], sep = "/")
  for (m in c("GMM", "MSGHD"))
    add(tag, "fit", m, r)
}

keys <- unique(keys)
seeds <- vapply(strsplit(keys, "|", fixed = TRUE), function(x)
  do.call(mk_seed, as.list(x)), integer(1))
dup <- duplicated(seeds) | duplicated(seeds, fromLast = TRUE)
if (any(dup))
  stop("semantic seed collision(s):\n", paste(keys[dup], seeds[dup], collapse = "\n"))
cat(length(keys), "semantic RNG keys audited; no collisions\n")
