# MixSIT_flip.R
#
# The axis-flip move. Part of the MixSIT estimator, not an add-on: fitMixSIT_ECM calls mixsit_flip_loop() before its final E-step, so every fit in the study --
# Sections 4, 5 and 6 alike -- uses the same estimator.

O_to_angles <- function(O, m) {
  r <- as.integer(det(O) < 0)
  P <- O; if (r == 1L) P[1, ] <- -P[1, ]
  IP <- diag(m) + P
  if (rcond(IP) < 1e-10) return(NULL)
  S <- (diag(m) - P) %*% solve(IP)
  list(angles = S[upper.tri(S)], reflection = r == 1L)
}

# Greedy: score every candidate flip, apply the best if it gains, repeat.

flip_polish_ref <- function(Y, params, lam_het, max_pass = 40, lam_U = 0, lam_O = 0,
                        mode = "U") {
  llp <- function(pp) {
    ex <- tryCatch(E_Step_cpp(Y, pp), error = function(e) NULL)
    if (is.null(ex)) return(-Inf)
    ex$loglik - compute_penalty_total(pp, lam_U, lam_O, 0, lam_het)
  }
  m <- ncol(Y); K <- length(params$pi)
  # re-solve one component's gamma on its current axes, with CM-2's own solver
  resolve_gamma <- function(th, z_k) {
    Dm <- th$O %*% backsolve(th$U, t(sweep(Y, 2, th$mu, "-")), transpose = TRUE)
    st <- log(pmax(th$gamma, 1.001) - 1)
    1 + exp(cm2_coord_cpp(Dm, z_k, th$nu, st,
                          lo = -20, hi = 4, ngrid = 49, nbisect = 60))
  }
  cur <- params; base <- llp(cur); n_applied <- 0L

  for (pass in seq_len(max_pass)) {
    ex <- tryCatch(E_Step_cpp(Y, cur), error = function(e) NULL)
    if (is.null(ex)) break
    best <- NULL; best_gain <- 1e-8

    if (mode == "full") {
      # one shared O: m candidate moves, each touching every component
      for (j in seq_len(m)) {
        Q <- cur
        Onew <- Q$components[[1]]$O; Onew[j, ] <- -Onew[j, ]
        ang <- O_to_angles(Onew, m); if (is.null(ang)) next
        for (k in seq_len(K)) {
          th <- Q$components[[k]]
          th$O <- Onew; th$angles <- ang$angles; th$reflection <- ang$reflection
          th$gamma <- resolve_gamma(th, ex$Z[, k])
          Q$components[[k]] <- th
        }
        g <- llp(Q) - base
        if (is.finite(g) && g > best_gain) { best_gain <- g; best <- Q }
      }
    } else {
      # VV and EV: O_k is private, so K * m candidate moves
      for (k in seq_len(K)) for (j in seq_len(m)) {
        Q <- cur; th <- Q$components[[k]]
        th$O[j, ] <- -th$O[j, ]
        ang <- O_to_angles(th$O, m); if (is.null(ang)) next
        th$angles <- ang$angles; th$reflection <- ang$reflection
        th$gamma <- resolve_gamma(th, ex$Z[, k])
        Q$components[[k]] <- th
        g <- llp(Q) - base
        if (is.finite(g) && g > best_gain) { best_gain <- g; best <- Q }
      }
    }

    if (is.null(best)) break
    cur <- best; base <- base + best_gain; n_applied <- n_applied + 1L
  }
  list(params = cur, n_flips = n_applied, llp = base)
}


.axis_terms <- function(D, nu, gm) {
  vapply(seq_along(nu), function(j) dskwet_vec_cpp(D[, j], nu[j], gm[j]),
         numeric(nrow(D)))
}
# latent coordinates, per-axis log-density terms, and the component's column of
# log(pi_k p_k(.)), all from one O(N m^2) pass
.comp_state <- function(Y, th, pi_k) {
  if (min(abs(diag(th$U))) < 1e-8) return(NULL)
  M <- tryCatch(backsolve(th$U, t(th$O)), error = function(e) NULL)
  if (is.null(M)) return(NULL)
  D <- sweep(Y, 2, th$mu, "-") %*% M
  Tm <- .axis_terms(D, th$nu, th$gamma)
  list(D = D, T = Tm,
       col = rowSums(Tm) - sum(log(diag(th$U))) + log(pi_k),
       cst = -sum(log(diag(th$U))) + log(pi_k))
}
.rowlse <- function(L) { mx <- do.call(pmax, as.data.frame(L))
  mx + log(rowSums(exp(L - mx))) }
# re-solve one axis only: cm2_coord_cpp on a single-row D
.solve_axis <- function(dj, z_k, nu_j, g_j) {
  th <- log(max(g_j, 1.001) - 1)
  1 + exp(cm2_coord_cpp(matrix(dj, nrow = 1), z_k, nu_j, th,
                        lo = -20, hi = 4, ngrid = 49, nbisect = 60))
}

flip_polish <- function(Y, params, lam_het, max_pass = 40, lam_U = 0, lam_O = 0,
                        mode = "U") {
  m <- ncol(Y); K <- length(params$pi)
  pen <- function(pp) compute_penalty_total(pp, lam_U, lam_O, 0, lam_het)
  cur <- params
  st <- lapply(seq_len(K), function(k) .comp_state(Y, cur$components[[k]], cur$pi[k]))
  if (any(vapply(st, is.null, TRUE)))
    return(flip_polish_ref(Y, params, lam_het, max_pass, lam_U, lam_O, mode))
  L <- vapply(st, `[[`, numeric(nrow(Y)), "col")
  if (!all(is.finite(.rowlse(L))))
    return(flip_polish_ref(Y, params, lam_het, max_pass, lam_U, lam_O, mode))
  base <- sum(.rowlse(L)) - pen(cur); n_applied <- 0L

  for (pass in seq_len(max_pass)) {
    Z <- exp(L - .rowlse(L))
    best <- NULL; best_gain <- 1e-8

    # gain from flipping axis j of component k, with only gamma_jk re-solved
    try_axis <- function(k, j) {
      dj <- -st[[k]]$D[, j]
      g  <- .solve_axis(dj, Z[, k], cur$components[[k]]$nu[j],
                        cur$components[[k]]$gamma[j])
      tj <- dskwet_vec_cpp(dj, cur$components[[k]]$nu[j], g)
      list(dj = dj, g = g, tj = tj,
           col = st[[k]]$col - st[[k]]$T[, j] + tj)
    }
    apply_axis <- function(Q, S, k, j, a, Onew, ang) {
      th <- Q$components[[k]]
      th$O <- Onew; th$angles <- ang$angles; th$reflection <- ang$reflection
      th$gamma[j] <- a$g
      Q$components[[k]] <- th
      S[[k]]$D[, j] <- a$dj; S[[k]]$T[, j] <- a$tj; S[[k]]$col <- a$col
      list(Q = Q, S = S)
    }

    cand <- if (mode == "full") seq_len(m) else NULL
    if (mode == "full") {
      for (j in cand) {
        Onew <- cur$components[[1]]$O; Onew[j, ] <- -Onew[j, ]
        ang <- O_to_angles(Onew, m); if (is.null(ang)) next
        Q <- cur; S <- st; Ln <- L
        for (k in seq_len(K)) {
          a <- try_axis(k, j)
          z <- apply_axis(Q, S, k, j, a, Onew, ang); Q <- z$Q; S <- z$S
          Ln[, k] <- a$col
        }
        v <- .rowlse(Ln)
        g <- if (all(is.finite(v))) sum(v) - pen(Q) - base else -Inf
        if (is.finite(g) && g > best_gain) { best_gain <- g
          best <- list(Q = Q, S = S, L = Ln) }
      }
    } else {
      for (k in seq_len(K)) for (j in seq_len(m)) {
        Onew <- cur$components[[k]]$O; Onew[j, ] <- -Onew[j, ]
        ang <- O_to_angles(Onew, m); if (is.null(ang)) next
        a <- try_axis(k, j)
        z <- apply_axis(cur, st, k, j, a, Onew, ang)
        Ln <- L; Ln[, k] <- a$col
        v <- .rowlse(Ln)
        g <- if (all(is.finite(v))) sum(v) - pen(z$Q) - base else -Inf
        if (is.finite(g) && g > best_gain) { best_gain <- g
          best <- list(Q = z$Q, S = z$S, L = Ln) }
      }
    }

    if (is.null(best)) break
    cur <- best$Q; st <- best$S; L <- best$L
    base <- base + best_gain; n_applied <- n_applied + 1L
  }
  list(params = cur, n_flips = n_applied, llp = base)
}

# Identifiability for the EE variant: one common permutation and one sign choice per shared axis.
canonicalise_shared_axes <- function(params) {
  K <- length(params$pi)
  m <- length(params$components[[1]]$gamma)
  O <- params$components[[1]]$O
  G <- vapply(params$components, function(cp) cp$gamma, numeric(m))
  Nu <- vapply(params$components, function(cp) cp$nu, numeric(m))
  LNu <- log(Nu)

  axis_score <- as.vector(LNu %*% params$pi)
  ord <- order(axis_score, method = "radix")
  changed_order <- !identical(ord, seq_len(m))
  G_new <- G
  Nu_new <- Nu
  if (changed_order) {
    O <- O[ord, , drop = FALSE]
    G_new <- G[ord, , drop = FALSE]
    Nu_new <- Nu[ord, , drop = FALSE]
  }

  need <- which(as.vector(log(G_new) %*% params$pi) < 0)
  if (length(need)) {
    O[need, ] <- -O[need, , drop = FALSE]
    G_new[need, ] <- 1 / G_new[need, , drop = FALSE]
  }
  if (!changed_order && !length(need)) return(params)
  ang <- O_to_angles(O, m)
  if (is.null(ang)) return(params)      # chart miss: leave the fit untouched
  for (k in seq_len(K)) {
    params$components[[k]]$gamma <- G_new[, k]
    params$components[[k]]$nu <- Nu_new[, k]
    params$components[[k]]$O <- O
    params$components[[k]]$angles <- ang$angles
    params$components[[k]]$reflection <- ang$reflection
  }
  params
}


# The gauge, applied at the end of an ECM iteration. A no-op for VV and EV, whose charts enforce their own domain by construction; for EE it is the projection
# described above. Separated from canonicalise_shared_axes() so the loops can call one function without testing the mode themselves.
gauge_canonicalise <- function(params, constraint_mode) {
  if (!identical(constraint_mode, "full")) return(params)
  out <- tryCatch(canonicalise_shared_axes(params), error = function(e) NULL)
  if (is.null(out)) params else out
}


# The estimator's flip stage: alternate flip-polish with ordinary ECM until no flip is accepted, or `outer` rounds are spent.
#
# Params in, params out. It deliberately does NOT run a final E-step or relabel the components -- fitMixSIT_ECM does both once, after this returns, so the
# labelling convention (order by max_j gamma_jk) is defined in exactly one place. Stage2/Section6_FlipLib.R wraps this for the paired base-vs-flip probes
# that need to polish an existing fit.
#
# `tol` is passed in rather than read from a config global so that the Stage 1 estimator has no Stage 2 dependency.
mixsit_flip_loop <- function(Y, P, lam_het, mode = "U", outer = 6, inner = 60,
                             lam_U = 0, lam_O = 0, tol = 1e-4) {
  if (identical(mode, "full")) return(list(params = P, n_flips = 0L))
  tot <- 0L
  for (it in seq_len(outer)) {
    fp <- tryCatch(flip_polish(Y, P, lam_het, lam_U = lam_U, lam_O = lam_O,
                               mode = mode), error = function(e) NULL)
    if (is.null(fp) || fp$n_flips == 0L) break
    tot <- tot + fp$n_flips
    P <- fp$params
    ex <- tryCatch(E_Step_cpp(Y, P), error = function(e) NULL)
    if (is.null(ex)) break
    prev <- ex$loglik - compute_penalty_total(P, lam_U, lam_O, 0, lam_het)
    for (t in seq_len(inner)) {
      st <- mixsit_aecm_step(Y, P, ex, lam_U, lam_O, 0, lam_het, mode)
      if (!isTRUE(st$accepted)) break
      P <- st$params
      ex <- st$expectations
      cur <- st$new_llp
      if (abs(cur - prev) / (abs(prev) + 1e-10) < tol) break
      prev <- cur
    }
  }
  list(params = P, n_flips = tot)
}
