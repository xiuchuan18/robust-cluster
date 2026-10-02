# Shared D3 recovery truth and exact one-to-one error matching.

SEP <- 6

# Exact square linear assignment by dynamic programming over column subsets. For the recovery design n = 10, so the O(n 2^n) cost is negligible and avoids
# both a package dependency and the finite-sample ambiguity of greedy matching. The returned integer vector maps each row to one distinct column.
solve_square_assignment <- function(cost) {
  cost <- as.matrix(cost)
  if (!is.numeric(cost) || nrow(cost) != ncol(cost) || !nrow(cost) ||
      any(!is.finite(cost)))
    stop("assignment cost must be a non-empty finite numeric square matrix")
  n <- nrow(cost)
  if (n >= 30L) stop("assignment dimension is too large for integer bit masks")
  bits <- bitwShiftL(1L, seq.int(0L, n - 1L))
  nstate <- bitwShiftL(1L, n)
  value <- rep(Inf, nstate)
  parent_mask <- parent_col <- rep(NA_integer_, nstate)
  value[1L] <- 0
  tol <- 1e-12

  for (mask in seq.int(0L, nstate - 2L)) {
    idx <- mask + 1L
    if (!is.finite(value[idx])) next
    row <- sum(bitwAnd(mask, bits) != 0L) + 1L
    for (col in seq_len(n)) if (bitwAnd(mask, bits[col]) == 0L) {
      new_mask <- bitwOr(mask, bits[col])
      candidate <- value[idx] + cost[row, col]
      new_idx <- new_mask + 1L
      if (candidate < value[new_idx] - tol) {
        value[new_idx] <- candidate
        parent_mask[new_idx] <- mask
        parent_col[new_idx] <- col
      }
    }
  }

  map <- integer(n)
  mask <- nstate - 1L
  for (row in n:1L) {
    idx <- mask + 1L
    col <- parent_col[idx]
    if (is.na(col)) stop("assignment backtracking failed")
    map[row] <- col
    mask <- parent_mask[idx]
  }
  map
}

# the truths
make_truth <- function(design, d, K, seed) {
  if (!identical(design, "D3")) stop("the release recovery design is D3")
  set.seed(seed)
  gsets <-
    list(seq(0.55, 1.75, length.out = d), seq(1.25, 2.25, length.out = d),
         seq(0.48, 1.44, length.out = d))
  nsets <- list(seq(4, 8, length.out = d), seq(3, 6, length.out = d),
                seq(5, 10, length.out = d))
  if (design == "D3") {          # permute (gamma, nu) jointly
    for (k in 1:K) { p <- sample.int(d)
      gsets[[k]] <- gsets[[k]][p]; nsets[[k]] <- nsets[[k]][p] }
  }
  U <- diag(d)
  if (design == "D3") {                # a real shared triangular U
    U <- diag(exp(seq(-0.35, 0.35, length.out = d)))
    U[upper.tri(U)] <- runif(d * (d - 1) / 2, -0.35, 0.35)
  }
  O <- lapply(1:K, function(k) qr.Q(qr(matrix(rnorm(d^2), d, d))))
  list(mu = list(rep(0, d), mu_at(d, 1:3, SEP), mu_at(d, 4:6, SEP)),
       nu = nsets, gm = gsets, U = U, O = O,
       lam = MIXSIT$lambda_hetero)
}

# error accounting
score <- function(P, cluster, lab, TR, d, K, design, N, r) {
  tab <- table(factor(lab, levels = 1:K), factor(cluster, levels = 1:K))
  # Rows are true components and columns are fitted components. Maximising the
  # total overlap is equivalent to minimising its negative.
  perm <- solve_square_assignment(-unclass(tab))
  out <- list()
  for (k in 1:K) {
    th <- P$components[[perm[k]]]
    Mh <- t(th$U) %*% t(th$O)                 # fitted observed-space axes
    Mt <- t(TR$U)  %*% t(TR$O[[k]])           # true ones
    Ch <- apply(Mh, 2, function(v) v / sqrt(sum(v^2)))
    Ct <- apply(Mt, 2, function(v) v / sqrt(sum(v^2)))
    g_true <- TR$gm[[k]]
    flip  <- ifelse(g_true >= 1, 1, -1)
    g_can <- ifelse(g_true >= 1, g_true, 1 / g_true)
    Ct_can <- sweep(Ct, 2, flip, "*")
    S_full <- t(Ct_can) %*% Ch
    S <- abs(S_full)
    S_sgn <- S_full
    # The projective angular loss is arccos(|q_true' q_fit|). Solve the whole one-to-one assignment rather than committing axes greedily.
    angle_cost <- acos(pmin(S, 1))
    map <- solve_square_assignment(angle_cost)
    me  <- vapply(1:d, function(j) fs_mean(g_true[j], TR$nu[[k]][j]), 0)
    loc <- TR$mu[[k]] - as.vector(Mt %*% me)
    # SIGN-AWARE gamma comparison. The axis match is made on |cos|, so a fitted axis may point along the canonicalised true axis or against it. Pointing
    # against it means the fitted axis represents the RECIPROCAL skew, so the target is 1/g_can, not g_can. Writing it this way is necessary now that EE
    # puts gamma on (0, inf) and its fitted values may legitimately fall below 1; for VV and EV, where gamma > 1 is enforced, it is also the stricter and more
    # honest comparison -- the previous version silently forgave an anti-aligned axis by comparing both sides on the gamma >= 1 representative.
    sgn <- ifelse(vapply(1:d, function(j) S_sgn[j, map[j]], 0) >= 0, 1, -1)
    out[[k]] <- data.frame(
      design = design, N = N, rep = r, k = k, axis = 1:d,
      matched_component = perm[k], matched_axis = map,
      component_matching = "global_maximum_overlap",
      axis_matching = "global_minimum_total_angle",
      mu_err  = sqrt(sum((th$mu - loc)^2)) / sqrt(d),
      lg_err  = log(th$gamma[map]) - sgn * log(g_can),
      lnu_err = log(th$nu[map]) - log(TR$nu[[k]]),
      g_can = g_can, g_hat = th$gamma[map],
      nu_true = TR$nu[[k]], nu_hat = th$nu[map],
      cos_abs = vapply(1:d, function(j) S[j, map[j]], 0),
      cos_sgn = vapply(1:d, function(j) sum(Ct_can[, j] * Ch[, map[j]]), 0),
      stringsAsFactors = FALSE)
  }
  R <- do.call(rbind, out)
  # U carries none of the axis symmetry, so it is compared directly
  R$U_err  <- sqrt(sum((P$components[[1]]$U - TR$U)^2)) / d
  R$pi_err <- max(abs(sort(P$pi) - 1 / K))
  R
}
