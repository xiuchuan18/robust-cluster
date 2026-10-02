# Regression tests for the direction-matched Wbar diagnostic.

setwd(normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "Stage2"),
  winslash = "/"))
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R")
}))

# Exact assignment must undo a signed permutation of observed axes.
d <- 5L
perm <- c(3L, 5L, 1L, 4L, 2L)
O <- diag(d)[perm, , drop = FALSE]
O[cbind(seq_len(d), perm)] <- O[cbind(seq_len(d), perm)] * c(-1, 1, -1, 1, 1)
component <- list(U = diag(d), O = O)
mt <- match_observed_axes(component)
stopifnot(identical(as.integer(mt$map), order(perm)),
          max(abs(mt$abs_cosine - 1)) < 1e-12)

# Posterior averaging must occur after component-specific matching.  Two
# components use different latent permutations but encode the same observed
# coordinate weights.
N <- 4L; K <- 2L
W <- array(NA_real_, c(N, d, K))
target <- outer(seq_len(N), 10 * seq_len(d), "+")
perm2 <- rev(seq_len(d))
O2 <- diag(d)[perm2, , drop = FALSE]
components <- list(list(U = diag(d), O = O), list(U = diag(d), O = O2))
maps <- lapply(components, match_observed_axes)
W[, maps[[1]]$map, 1] <- target
W[, maps[[2]]$map, 2] <- target
Z <- cbind(c(.2, .4, .6, .8), c(.8, .6, .4, .2))
out <- matrix(0, N, d)
for (k in seq_len(K)) out <- out + Z[, k] * W[, maps[[k]]$map, k]
stopifnot(max(abs(out - target)) < 1e-12)

cat("Wbar axis-matching tests passed\n")
