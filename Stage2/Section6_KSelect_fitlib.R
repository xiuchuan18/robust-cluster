# Section6_KSelect_fitlib.R
#
# fit_at_K(): one method, one K, returning the raw log-likelihood, the posterior and the partition. Shared by the main K-selection grid and the difficulty
# ladder so that the two cannot drift apart -- a divergence in how a method is fitted between them would make the ladder uninterpretable as an extension of the grid.
#
# The MixSIT branch uses the axis-flip move, because fitMixSIT_ECM applies it by default and the estimator must be the same one Sections 4 and 5 report.
#
# Requires Section4_Lib.R (for .quiet and the MixSIT fitter) and MixGHD.
#
# `tag` and `rep` seed each method's own stream; see fit_multistart.
fit_at_K <- function(method, Y, K, d, tag = NULL, rep = NULL) {
  n <- nrow(Y)
  code <- if (method %in% c("GMM", "MSGHD")) method else method
  r <- tryCatch(fit_multistart(code, Y, K, R_RESTART, COMPETITOR_MAX_ITER,
                               tag = tag, rep = if (is.null(rep)) NULL else
                                 paste(rep, "K", K, sep = "_")),
                error = function(e) NULL)
  if (is.null(r) || is.null(r$fit)) return(NULL)
  f <- r$fit

  if (method == "MSGHD") {
    ll <- tryCatch(as.numeric(tail(f@loglik, 1)), error = function(e) NA_real_)
    if (!is.finite(ll)) return(NULL)
    return(list(ll = ll, z = model_posterior(f, method, n, K),
                cl = r$cluster))
  }
  if (method == "GMM")
    return(list(ll = as.numeric(f$loglik), z = f$z, cl = r$cluster))

  # MixSIT: the RAW log-likelihood, not the penalised objective the algorithm ascends, so the criteria mean the same thing as they do for the competitors.
  list(ll = as.numeric(f$loglik_obs), z = f$posterior, cl = r$cluster)
}
