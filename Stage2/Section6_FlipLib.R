# Section6_FlipLib.R
#
# The axis-flip move now lives in Stage1/MixSIT_flip.R and is applied by
# fitMixSIT_ECM itself (argument `flips`, on by default), so that Sections 4, 5
# and 6 all report the same estimator. Nothing is defined here any more; the
# file remains only as the two wrappers the paired probes are written against.

polish_from <- function(Y, P, K, lam_het, mode = "U", outer = 6, inner = 60,
                        lam_U = 0, lam_O = 0) {
  fp <- mixsit_flip_loop(Y, P, lam_het, mode = mode, outer = outer,
                         inner = inner, lam_U = lam_U, lam_O = lam_O,
                         tol = MIXSIT$tol)
  P <- fp$params
  ex <- E_Step_cpp(Y, P)
  if (length(P$pi) > 1) {
    o <- order(sapply(P$components, function(c) max(c$gamma)))
    P$components <- P$components[o]; P$pi <- P$pi[o]
    ex$Z <- ex$Z[, o, drop = FALSE]
  }
  list(params = P, cluster = apply(ex$Z, 1, which.max), n_flips = fp$n_flips,
       loglik_obs = as.numeric(ex$loglik),
       llp = ex$loglik - compute_penalty_total(P, lam_U, lam_O, 0, lam_het))
}

# A single-start fit with the flip stage, at the study's restart budget. Kept because the Section 6 probes call it by name; it is now exactly what
# fitMixSIT_ECM does by default.
fit_with_flips <- function(Y, K, lam_het, mode = "U", outer = 6, inner = 60,
                           lam_U = 0, lam_O = 0) {
  fitMixSIT_ECM(Y, K = K, constraint_mode = mode, n_short = 1, short_iter = 1,
                init_methods = "kmeans", max_iter = 500, tol = MIXSIT$tol,
                lambda_U = lam_U, lambda_O = lam_O, lambda_magnitude = 0,
                lambda_hetero = lam_het,
                flips = TRUE, flip_outer = outer, flip_inner = inner)
}
