args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
lib <- file.path(root, ".r-lib")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

source(file.path(root, "Stage1", "MixSIT_rcpp.R"))
source(file.path(root, "Stage1", "MixSIT_flip.R"))

same_matrix <- function(a, b, tol = 1e-10)
  isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))

set.seed(290926)
n <- 80L; m <- 4L
Y <- rbind(
  matrix(stats::rt(n * m, df = 4), n, m),
  sweep(matrix(stats::rt(n * m, df = 6), n, m), 2,
        c(3, -2, 1, 0), "+")
)
bad <- matrix(stats::runif(length(Y)) < 0.04, nrow(Y), ncol(Y))
Y[bad] <- Y[bad] + stats::rnorm(sum(bad), 0, 6)

# EV and EE must be feasible before their first E-step.
for (mode in c("U", "full")) {
  p <- initialize_mixsit(Y, 2, "kmeans", mode)
  stopifnot(same_matrix(p$components[[1]]$U, p$components[[2]]$U))
  stopifnot(is.finite(E_Step_cpp(Y, p)$loglik))
}

# EE uses one common permutation.  The canonicalization must preserve the
# likelihood while allowing a component-specific tail ranking.
p <- initialize_mixsit(Y, 2, "kmeans", "full")
p$components[[1]]$nu <- c(20, 3, 8, 11)
p$components[[2]]$nu <- c(2, 15, 6, 9)
p$components[[1]]$gamma <- c(0.7, 2.0, 0.8, 1.1)
p$components[[2]]$gamma <- c(0.6, 0.9, 1.4, 0.75)
ang <- c(0.17, -0.23, 0.11, 0.19, -0.07, 0.13)
O0 <- build_O_R(ang, m, FALSE)
for (k in 1:2) {
  p$components[[k]]$O <- O0
  p$components[[k]]$angles <- ang
}
ll0 <- E_Step_cpp(Y, p)$loglik
q <- canonicalise_shared_axes(p)
ll1 <- E_Step_cpp(Y, q)$loglik
stopifnot(abs(ll1 - ll0) < 1e-8)
scores <- as.vector(log(vapply(q$components, function(x) x$nu, numeric(m))) %*% q$pi)
stopifnot(all(diff(scores) >= -1e-12))
stopifnot(any(diff(q$components[[1]]$nu) < 0) ||
          any(diff(q$components[[2]]$nu) < 0))
stopifnot(all(vapply(q$components, function(x)
  all(is.finite(x$nu)) && all(x$nu > 0), logical(1))))

# An exact M1 tie is stable: neither higher tail moments nor gamma magnitudes
# are used to reorder the tied axes.
pt <- p
pt$pi <- c(0.5, 0.5)
for (k in 1:2) {
  pt$components[[k]]$nu <- c(3, 3, 8, 12)
  pt$components[[k]]$gamma <- c(3, 1.1, 1.4, 1.6)
  pt$components[[k]]$O <- O0
  pt$components[[k]]$angles <- ang
}
llt0 <- E_Step_cpp(Y, pt)$loglik
qt <- canonicalise_shared_axes(pt)
stopifnot(abs(E_Step_cpp(Y, qt)$loglik - llt0) < 1e-8)
stopifnot(qt$components[[1]]$gamma[1] > qt$components[[1]]$gamma[2])
tie_scores <- as.vector(log(vapply(qt$components,
                                   function(x) x$nu, numeric(m))) %*% qt$pi)
stopifnot(tie_scores[1] == tie_scores[2])

# The shared-axis gauge also applies at K=1.
p1 <- initialize_mixsit(Y[1:n, , drop = FALSE], 1, "kmeans", "full")
p1$components[[1]]$gamma <- c(0.5, 2, 0.8, 1.1)
p1$components[[1]]$O <- O0
p1$components[[1]]$angles <- ang
q1 <- gauge_canonicalise(p1, "full")
stopifnot(all(log(q1$components[[1]]$gamma) >= -1e-12))
stopifnot(abs(E_Step_cpp(Y[1:n, , drop = FALSE], p1)$loglik -
              E_Step_cpp(Y[1:n, , drop = FALSE], q1)$loglik) < 1e-8)

# A short difficult fit in every structural mode must be monotone at the
# observed penalized objective and must not need the whole-iteration rejection.
for (mode in c("none", "U", "full")) {
  set.seed(8100 + match(mode, c("none", "U", "full")))
  fit <- fitMixSIT_ECM(
    Y, 2, constraint_mode = mode,
    n_short = 1, short_iter = 1, init_methods = "kmeans",
    max_iter = 15, tol = 1e-8, flips = FALSE)
  stopifnot(all(is.finite(fit$loglik)))
  stopifnot(all(diff(fit$loglik) >= -1e-7))
  stopifnot(identical(fit$ascent_rejections, 0L))
  if (mode %in% c("U", "full"))
    stopifnot(same_matrix(fit$params$components[[1]]$U,
                          fit$params$components[[2]]$U))
  if (mode == "full") {
    stopifnot(same_matrix(fit$params$components[[1]]$O,
                          fit$params$components[[2]]$O))
    sc <- as.vector(log(vapply(fit$params$components,
                              function(x) x$nu, numeric(m))) %*% fit$params$pi)
    stopifnot(all(diff(sc) >= -1e-10))
  }
}

# Section 7's penalty-sensitivity library must preserve the two-cycle AECM
# interface.  In particular, it may replace CM_Step_Nu but must not recreate the
# retired combined shape step that reused stale latent-scale moments.
source(file.path(root, "Stage2", "Section6_PenaltyLib.R"))
install_penalty("log", 2)
stopifnot(exists("CM_Step_Nu", inherits = FALSE),
          !exists("CM_Step_Shape", inherits = FALSE),
          check_penalty_gradient())
for (mode in c("U", "full")) {
  set.seed(620L + match(mode, c("U", "full")))
  fit <- fitMixSIT_ECM(
    Y, K = 2, constraint_mode = mode,
    n_short = 1, short_iter = 1, init_methods = "kmeans",
    max_iter = 5, tol = 1e-8, lambda_U = 0, lambda_O = 0,
    lambda_magnitude = 0, lambda_hetero = 2, flips = FALSE)
  stopifnot(all(is.finite(fit$loglik)),
            all(diff(fit$loglik) >= -1e-6),
            identical(fit$ascent_rejections, 0L))
}
cat("penalty-library compatibility tests passed\n")
cat("final method revision tests passed\n")
