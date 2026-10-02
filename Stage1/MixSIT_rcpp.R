# MixSIT_rcpp.R
#
# R-side ECM driver for the MixSIT mixture model.


suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(mclust)
})

# Load the precompiled C++ backend from the MixSITcore package. (Rcpp::sourceCpp() does not work reliably inside parallel workers on
# Windows; the canonical fix is to install SIT_core.cpp as a real R package. Run Build_MixSITcore_Package.R once to create it.)
if (!requireNamespace("MixSITcore", quietly = TRUE)) {
  stop(
    "Package 'MixSITcore' is not installed. Run Build_MixSITcore_Package.R ",
    "from the Stage1 directory first."
  )
}
suppressPackageStartupMessages(library(MixSITcore))

# Helper: penalty total (P_nu + P_O + P_U) summed over components.
compute_penalty_total <- function(params,
                                  lambda_U, lambda_O,
                                  lambda_magnitude, lambda_hetero) {
  total <- 0.0
  K <- length(params$pi)
  m <- length(params$components[[1]]$nu)
  for (k in seq_len(K)) {
    comp <- params$components[[k]]

    # P_nu
    log_nu <- log(comp$nu)
    p_mag <- if (lambda_magnitude > 0)
      lambda_magnitude * sum((log_nu - log(10))^2) else 0
    p_het <- if (lambda_hetero > 0)
      lambda_hetero * sum((log_nu - mean(log_nu))^2) else 0

    # P_O (uses angles vector; matches the C++ optimizer's penalty)
    angles <- if (!is.null(comp$angles)) comp$angles else numeric(0)
    p_O <- if (length(angles) > 0) lambda_O * sum(angles^2) else 0

    # P_U (strict upper triangle)
    U <- comp$U
    p_U <- if (m > 1) lambda_U * sum(U[upper.tri(U)]^2) else 0

    total <- total + p_mag + p_het + p_O + p_U
  }
  total
}

# CM-Step Pi: closed-form update of mixing proportions.
CM_Step_Pi <- function(Z_hat) {
  N <- nrow(Z_hat)
  pi_new <- colSums(Z_hat) / N
  pi_new[pi_new < 1e-10] <- 1e-10
  pi_new / sum(pi_new)
}

# CM-Step Shape: update gamma_k (AECM) then nu_k.

CM2_THETA_RANGE <- c(-20, 4)
CM2_NGRID       <- 49
CM2_NBISECT     <- 60

# Guarded analytic-gradient solvers for the CM-2 and CM-3 subproblems.

.cm2_vg <- function(D_mat, z_k, nu) {
  cache_par <- NULL; cache_res <- NULL
  ev <- function(th) {
    if (is.null(cache_par) || !identical(th, cache_par)) {
      g <- 1.0 + exp(th)
      r <- obj_grad_gamma_cpp(D_mat, z_k, nu, g)
      # chain rule for gamma = 1 + exp(theta): d gamma / d theta = gamma - 1
      cache_res <<- list(value = r$value, grad = r$grad * (g - 1.0))
      cache_par <<- th
    }
    cache_res
  }
  list(fn = function(th) ev(th)$value, gr = function(th) ev(th)$grad)
}

# Value and gradient of the CM-3 objective in the increment parameterisation.
.cm3_vg <- function(sum_z, tilde_d, lambda_magnitude, lambda_hetero) {
  m <- length(tilde_d)
  cache_par <- NULL; cache_res <- NULL
  ev <- function(dl) {
    if (is.null(cache_par) || !identical(dl, cache_par)) {
      nu <- cumsum(exp(dl)); v <- nu / 2; lnu <- log(nu)
      val <- -(sum_z * sum(v * log(v) - lgamma(v)) + sum(v * tilde_d))
      gnu <- -(sum_z / 2) * (log(v) + 1 - digamma(v)) - tilde_d / 2
      if (lambda_magnitude > 0) {
        val <- val + lambda_magnitude * sum((lnu - log(10))^2)
        gnu <- gnu + 2 * lambda_magnitude * (lnu - log(10)) / nu
      }
      if (lambda_hetero > 0) {
        lbar <- mean(lnu)
        val <- val + lambda_hetero * sum((lnu - lbar)^2)
        gnu <- gnu + 2 * lambda_hetero * (lnu - lbar) / nu
      }
      # d/d delta_t = exp(delta_t) * sum_{j >= t} d/d nu_j
      gdl <- exp(dl) * rev(cumsum(rev(gnu)))
      cache_res <<- list(value = val, grad = gdl)
      cache_par <<- dl
    }
    cache_res
  }
  list(fn = function(dl) ev(dl)$value, gr = function(dl) ev(dl)$grad)
}

# EE has one common latent-axis permutation, not one independent permutation per component.  Its tail indices are therefore positive but otherwise free;
# a single common axis ordering is imposed later by canonicalise_shared_axes().
.cm3_free_vg <- function(sum_z, tilde_d, lambda_magnitude, lambda_hetero) {
  cache_par <- NULL; cache_res <- NULL
  ev <- function(eta) {
    if (is.null(cache_par) || !identical(eta, cache_par)) {
      nu <- exp(eta); v <- nu / 2; lnu <- eta
      val <- -(sum_z * sum(v * log(v) - lgamma(v)) + sum(v * tilde_d))
      gnu <- -(sum_z / 2) * (log(v) + 1 - digamma(v)) - tilde_d / 2
      if (lambda_magnitude > 0) {
        val <- val + lambda_magnitude * sum((lnu - log(10))^2)
        gnu <- gnu + 2 * lambda_magnitude * (lnu - log(10)) / nu
      }
      if (lambda_hetero > 0) {
        lbar <- mean(lnu)
        val <- val + lambda_hetero * sum((lnu - lbar)^2)
        gnu <- gnu + 2 * lambda_hetero * (lnu - lbar) / nu
      }
      cache_res <<- list(value = val, grad = gnu * nu)
      cache_par <<- eta
    }
    cache_res
  }
  list(fn = function(eta) ev(eta)$value, gr = function(eta) ev(eta)$grad)
}

.cm2_solve <- function(par0, vg, fd_fn, control) {
  j0 <- vg$fn(par0)
  o <- tryCatch(nlminb(par0, vg$fn, vg$gr, control = control),
                error = function(e) NULL)
  if (!is.null(o) && is.finite(o$objective) && o$objective <= j0) return(o)
  o2 <- tryCatch(nlminb(par0, fd_fn, control = control), error = function(e) NULL)
  if (!is.null(o2) && is.finite(o2$objective) && o2$objective <= j0) return(o2)
  list(par = par0, objective = j0)
}

# CM-2 for the EE variant: gamma free on (0, infinity).
#

cm2_gamma_free <- function(D_mat, z_k, nu, gamma0) {
  m <- nrow(D_mat)
  solve_side <- function(D, g0) {
    t0 <- log(pmax(g0, 1.001) - 1)
    1 + exp(cm2_coord_cpp(D, z_k, nu, t0,
                          lo = CM2_THETA_RANGE[1], hi = CM2_THETA_RANGE[2],
                          ngrid = CM2_NGRID, nbisect = CM2_NBISECT))
  }
  g_up <- solve_side( D_mat, gamma0)          # the gamma >= 1 branch
  g_dn <- solve_side(-D_mat, 1 / gamma0)      # the gamma <= 1 branch, as 1/gamma
  # per-axis objective, so each axis takes its own better orientation
  fj <- function(dj, nu_j, g) -sum(z_k * dskwet_vec_cpp(dj, nu_j, g))
  vapply(seq_len(m), function(j) {
    dj <- D_mat[j, ]
    if (fj(dj, nu[j], 1 / g_dn[j]) < fj(dj, nu[j], g_up[j])) 1 / g_dn[j] else g_up[j]
  }, numeric(1))
}

CM_Step_Gamma <- function(Y, params, expectations, constraint_mode = "U") {
  K <- length(params$pi)
  m <- ncol(Y)
  Z_hat <- expectations$Z
  new_components <- params$components

  for (k in seq_len(K)) {
    theta <- new_components[[k]]
    z_k <- Z_hat[, k]
    if (sum(z_k) < 1e-6) next

    Y_centered <- sweep(Y, 2, theta$mu, "-")
    Temp <- backsolve(theta$U, t(Y_centered), transpose = TRUE)
    D_mat <- theta$O %*% Temp
    safe_gamma <- pmax(theta$gamma, 1.001)
    start_theta <- log(safe_gamma - 1.0)

    # Cycle 1 of the AECM integrates out W.  The objective is separable over axes and retains the incoming coordinate as an exact non-worsening guard.
    new_components[[k]]$gamma <- if (identical(constraint_mode, "full")) {
      cm2_gamma_free(D_mat, z_k, theta$nu, theta$gamma)
    } else {
      theta_new <- cm2_coord_cpp(D_mat, z_k, theta$nu, start_theta,
                                 lo = CM2_THETA_RANGE[1], hi = CM2_THETA_RANGE[2],
                                 ngrid = CM2_NGRID, nbisect = CM2_NBISECT)
      1.0 + exp(theta_new)
    }
  }
  new_components
}

CM_Step_Nu <- function(params, expectations,
                       lambda_magnitude, lambda_hetero,
                       constraint_mode = "U") {
  K <- length(params$pi)
  m <- length(params$components[[1]]$nu)
  Z_hat <- expectations$Z
  W_hat <- expectations$W
  Kappa_hat <- expectations$Kappa
  new_components <- params$components

  for (k in seq_len(K)) {
    theta <- new_components[[k]]
    z_k <- Z_hat[, k]
    sum_z <- sum(z_k)
    if (sum_z < 1e-6) next

    tilde_d <- vapply(seq_len(m), function(j)
      sum(z_k * (Kappa_hat[, j, k] - W_hat[, j, k])), numeric(1))

    if (identical(constraint_mode, "full")) {
      # EE has only a common axis permutation.  Do not sort each component's tail indices independently: positivity is the only within-component
      # constraint, and one common permutation is fixed after the structural CM.
      obj_nu_free <- function(eta) {
        nu_vec <- exp(eta); v_half <- nu_vec / 2
        neg_obj <- -(sum_z * sum(v_half * log(v_half) - lgamma(v_half)) +
                     sum(v_half * tilde_d))
        if (lambda_magnitude > 0)
          neg_obj <- neg_obj + lambda_magnitude * sum((eta - log(10))^2)
        if (lambda_hetero > 0)
          neg_obj <- neg_obj + lambda_hetero * sum((eta - mean(eta))^2)
        neg_obj
      }
      eta0 <- log(pmax(theta$nu, 1e-8))
      opt_nu <- .cm2_solve(
        eta0,
        .cm3_free_vg(sum_z, tilde_d, lambda_magnitude, lambda_hetero),
        obj_nu_free,
        list(iter.max = 100, eval.max = 100, rel.tol = 1e-6))
      new_components[[k]]$nu <- exp(opt_nu$par)
    } else {
      obj_nu_inc <- function(delta_vec) {
        nu_vec <- cumsum(exp(delta_vec)); v_half <- nu_vec / 2
        neg_obj <- -(sum_z * sum(v_half * log(v_half) - lgamma(v_half)) +
                     sum(v_half * tilde_d))
        log_nu_vec <- log(nu_vec)
        if (lambda_magnitude > 0)
          neg_obj <- neg_obj + lambda_magnitude * sum((log_nu_vec - log(10))^2)
        if (lambda_hetero > 0)
          neg_obj <- neg_obj + lambda_hetero *
            sum((log_nu_vec - mean(log_nu_vec))^2)
        neg_obj
      }
      sorted_nu <- sort(theta$nu)
      delta_init <- numeric(m)
      delta_init[1] <- log(max(sorted_nu[1], 2))
      if (m > 1)
        delta_init[2:m] <- log(pmax(diff(sorted_nu), 0.1))
      opt_nu <- .cm2_solve(
        delta_init,
        .cm3_vg(sum_z, tilde_d, lambda_magnitude, lambda_hetero),
        obj_nu_inc,
        list(iter.max = 100, eval.max = 100, rel.tol = 1e-6))
      new_components[[k]]$nu <- cumsum(exp(opt_nu$par))
    }

    # The solver guards its encoded start, which may differ from theta$nu after initialization floors. Guard the actual incoming parameter too.
    # Evaluate both values in nu space, including both penalties, so a failed or worsening candidate retains the original nu exactly (also for EE).
    objective_actual_nu <- function(nu_vec) {
      if (any(!is.finite(nu_vec)) || any(nu_vec <= 0)) return(Inf)
      v_half <- nu_vec / 2
      log_nu_vec <- log(nu_vec)
      value <- -(sum_z * sum(v_half * log(v_half) - lgamma(v_half)) +
                 sum(v_half * tilde_d))
      if (lambda_magnitude > 0)
        value <- value + lambda_magnitude * sum((log_nu_vec - log(10))^2)
      if (lambda_hetero > 0)
        value <- value + lambda_hetero *
          sum((log_nu_vec - mean(log_nu_vec))^2)
      value
    }
    old_objective <- objective_actual_nu(theta$nu)
    new_objective <- objective_actual_nu(new_components[[k]]$nu)
    if (!is.finite(old_objective) || !is.finite(new_objective) ||
        new_objective > old_objective)
      new_components[[k]]$nu <- theta$nu
  }
  new_components
}

# Guarded analytic-gradient solver for the CM-4 subproblems.
.cm4_vg <- function(f) {
  cache_par <- NULL; cache_res <- NULL
  ev <- function(p) {
    if (is.null(cache_par) || !identical(p, cache_par)) {
      cache_res <<- f(p); cache_par <<- p
    }
    cache_res
  }
  list(fn = function(p) ev(p)$value, gr = function(p) ev(p)$grad)
}


.cm4_solve <- function(par0, vg, fd_fn, control) {
  j0 <- vg$fn(par0)
  o <- tryCatch({
    r <- optim(par0, vg$fn, vg$gr, method = "L-BFGS-B",
               control = list(maxit = control$iter.max, factr = 1e7))
    list(par = r$par, objective = r$value)
  }, error = function(e) NULL)
  if (!is.null(o) && is.finite(o$objective) && o$objective <= j0) return(o)
  # the trust-region method is the first fallback, the differenced one the last
  o2 <- tryCatch(nlminb(par0, vg$fn, vg$gr, control = control),
                 error = function(e) NULL)
  if (!is.null(o2) && is.finite(o2$objective) && o2$objective <= j0) return(o2)
  o3 <- tryCatch(nlminb(par0, fd_fn, control = control), error = function(e) NULL)
  if (!is.null(o3) && is.finite(o3$objective) && o3$objective <= j0) return(o3)
  list(par = par0, objective = j0)
}

# CM-Step Space: update (mu_k, U_k, S_k) under the chosen constraint_mode. Three modes: "none" (VV, per-component), "U" (EV, shared U with per-k O),
# "full" (EE, shared U and shared O). Dual-path (rotation/reflection) search.
build_U_R <- function(u_pars, m) {
  u_diag <- exp(pmin(pmax(u_pars[1:m], -10), 10))
  U <- diag(u_diag, m, m)
  if (m > 1) {
    idx <- m
    for (j in seq_len(m)) {
      for (i in seq_len(j - 1)) {
        if (idx + 1 <= length(u_pars)) U[i, j] <- u_pars[idx + 1]
        idx <- idx + 1
      }
    }
  }
  U
}

build_O_R <- function(ang_pars, m, use_reflection = FALSE) {
  vector_to_orthogonal_cpp(ang_pars, m, use_reflection)
}

CM_Step_Space <- function(Y, params, expectations,
                          lambda_U, lambda_O, constraint_mode) {
  K <- length(params$pi)
  m <- ncol(Y)
  new_components <- params$components

  len_angles <- if (m > 1) m * (m - 1) / 2 else 0
  len_U_diag <- m
  len_U_vec  <- if (m > 1) m * (m - 1) / 2 else 0

  opt_control <- list(iter.max = 100, eval.max = 200, rel.tol = 1e-6)

  if (constraint_mode == "none") {
    # Per-component (mu_k, U_k, O_k); dual-path on reflection.
    for (k in seq_len(K)) {
      theta <- new_components[[k]]
      curr_ang <- if (!is.null(theta$angles)) theta$angles else rep(0, len_angles)
      par_init <- c(theta$mu,
                    log(diag(theta$U)),
                    if (m > 1) theta$U[upper.tri(theta$U)] else numeric(0),
                    curr_ang)

      # PERF: slice the posterior and the latent-weight cube ONCE per CM-4 call.
      # These used to sit inside the closures, so R re-created an N x m copy of
      # W[, , k] on every one of the ~1000-5000 objective evaluations nlminb makes.
      z_k <- expectations$Z[, k]
      w_k <- expectations$W[, , k]
      opt_rot <- .cm4_solve(par_init,
        .cm4_vg(function(p) obj_grad_hetero_cpp(p, Y, theta$gamma, z_k, w_k,
                                                lambda_U, lambda_O, FALSE)),
        function(p) obj_space_hetero_cpp(p, Y, theta, z_k, w_k,
                                         lambda_U, lambda_O, FALSE), opt_control)
      opt_ref <- .cm4_solve(par_init,
        .cm4_vg(function(p) obj_grad_hetero_cpp(p, Y, theta$gamma, z_k, w_k,
                                                lambda_U, lambda_O, TRUE)),
        function(p) obj_space_hetero_cpp(p, Y, theta, z_k, w_k,
                                         lambda_U, lambda_O, TRUE), opt_control)

      if (opt_ref$objective < opt_rot$objective) {
        best_par <- opt_ref$par;  use_ref <- TRUE
      } else {
        best_par <- opt_rot$par;  use_ref <- FALSE
      }

      idx <- 1
      new_components[[k]]$mu <- best_par[idx:(idx + m - 1)]; idx <- idx + m
      u_pars <- best_par[idx:(idx + len_U_diag + len_U_vec - 1)]
      idx <- idx + len_U_diag + len_U_vec
      new_components[[k]]$U <- build_U_R(u_pars, m)
      ang <- if (m > 1) best_par[idx:(idx + len_angles - 1)] else numeric(0)
      new_components[[k]]$angles <- ang
      new_components[[k]]$reflection <- use_ref
      new_components[[k]]$O <- build_O_R(ang, m, use_ref)
    }

  } else if (constraint_mode == "full") {
    # Shared U and shared O; dual-path on reflection of the shared O.
    u0 <- new_components[[1]]$U
    par_init <- numeric(0)
    for (k in seq_len(K)) par_init <- c(par_init, new_components[[k]]$mu)
    par_init <- c(par_init, log(diag(u0)),
                  if (m > 1) u0[upper.tri(u0)] else numeric(0))
    ang0 <- if (!is.null(new_components[[1]]$angles))
      new_components[[1]]$angles else rep(0, len_angles)
    par_init <- c(par_init, ang0)

    Zh <- expectations$Z; Wh <- expectations$W        # PERF: hoisted
    opt_rot <- .cm4_solve(par_init,
      .cm4_vg(function(p) obj_grad_full_cpp(p, Y, new_components, Zh, Wh,
                                            lambda_U, lambda_O, FALSE)),
      function(p) obj_space_full_cpp(p, Y, new_components, Zh, Wh,
                                     lambda_U, lambda_O, FALSE), opt_control)
    opt_ref <- .cm4_solve(par_init,
      .cm4_vg(function(p) obj_grad_full_cpp(p, Y, new_components, Zh, Wh,
                                            lambda_U, lambda_O, TRUE)),
      function(p) obj_space_full_cpp(p, Y, new_components, Zh, Wh,
                                     lambda_U, lambda_O, TRUE), opt_control)

    if (opt_ref$objective < opt_rot$objective) {
      best_par <- opt_ref$par; use_ref <- TRUE
    } else {
      best_par <- opt_rot$par; use_ref <- FALSE
    }

    idx <- 1
    for (k in seq_len(K)) {
      new_components[[k]]$mu <- best_par[idx:(idx + m - 1)]; idx <- idx + m
    }
    u_pars <- best_par[idx:(idx + len_U_diag + len_U_vec - 1)]
    idx <- idx + len_U_diag + len_U_vec
    U_shared <- build_U_R(u_pars, m)
    ang <- if (m > 1) best_par[idx:(idx + len_angles - 1)] else numeric(0)
    O_shared <- build_O_R(ang, m, use_ref)
    for (k in seq_len(K)) {
      new_components[[k]]$U <- U_shared
      new_components[[k]]$O <- O_shared
      new_components[[k]]$angles <- ang
      new_components[[k]]$reflection <- use_ref
    }

  } else if (constraint_mode == "U") {
    # Shared U with per-component O_k.
    # Stage 1: pre-scan to determine reflection direction per cluster (with current U fixed).
    # Stage 2: joint optimization over shared U and all (mu_k, O_k).
    current_U <- new_components[[1]]$U
    ref_vec <- logical(K)
    U_diag_log <- log(diag(current_U))

    for (k in seq_len(K)) {
      theta <- new_components[[k]]
      curr_ang <- if (!is.null(theta$angles)) theta$angles else rep(0, len_angles)

      z_k <- expectations$Z[, k]
      w_k <- expectations$W[, , k]
      u_upper <- if (m > 1) current_U[upper.tri(current_U)] else numeric(0)
      nA  <- len_angles
      gsl <- if (nA > 0) (2 * m + nA + 1):(2 * m + 2 * nA) else integer(0)
      mk_vg <- function(is_ref) {
        cp <- NULL; cr <- NULL
        ev <- function(x) {
          if (is.null(cp) || !identical(x, cp)) {
            p_vec <- c(theta$mu, U_diag_log, u_upper, x)
            r <- obj_grad_hetero_cpp(p_vec, Y, theta$gamma, z_k, w_k,
                                     lambda_U, lambda_O, is_ref)
            cr <<- list(value = r$value, grad = r$grad[gsl]); cp <<- x
          }
          cr
        }
        list(fn = function(x) ev(x)$value, gr = function(x) ev(x)$grad)
      }
      vg_rot <- mk_vg(FALSE); vg_ref <- mk_vg(TRUE)
      res_rot <- nlminb(start = curr_ang, objective = vg_rot$fn,
                        gradient = vg_rot$gr,
                        control = list(iter.max = 20, eval.max = 50))
      res_ref <- nlminb(start = curr_ang, objective = vg_ref$fn,
                         gradient = vg_ref$gr,
                         control = list(iter.max = 20, eval.max = 50))

      # Include the incoming chart and coordinates explicitly.  The pre-scan chooses a disconnected chart for the joint solve that follows, so it
      # needs the same non-worsening protection as every other structural step.
      curr_ref <- isTRUE(theta$reflection)
      best_val <- if (curr_ref) vg_ref$fn(curr_ang) else vg_rot$fn(curr_ang)
      best_ang <- curr_ang
      best_ref <- curr_ref
      if (is.finite(res_rot$objective) && res_rot$objective <= best_val) {
        best_val <- res_rot$objective; best_ang <- res_rot$par; best_ref <- FALSE
      }
      if (is.finite(res_ref$objective) && res_ref$objective <= best_val) {
        best_val <- res_ref$objective; best_ang <- res_ref$par; best_ref <- TRUE
      }
      ref_vec[k] <- best_ref
      new_components[[k]]$angles <- best_ang
    }

    par_init <- numeric(0)
    for (k in seq_len(K)) par_init <- c(par_init, new_components[[k]]$mu)
    for (k in seq_len(K)) par_init <- c(par_init, new_components[[k]]$angles)
    par_init <- c(par_init, log(diag(current_U)),
                  if (m > 1) current_U[upper.tri(current_U)] else numeric(0))

    Zh <- expectations$Z; Wh <- expectations$W        # PERF: hoisted
    opt_global <- .cm4_solve(par_init,
      .cm4_vg(function(p) obj_grad_u_only_cpp(p, Y, new_components, Zh, Wh,
                                              lambda_U, lambda_O, ref_vec)),
      function(p) obj_space_u_only_cpp(p, Y, new_components, Zh, Wh,
                                       lambda_U, lambda_O, ref_vec), opt_control)

    idx <- 1
    for (k in seq_len(K)) {
      new_components[[k]]$mu <- opt_global$par[idx:(idx + m - 1)]; idx <- idx + m
    }
    for (k in seq_len(K)) {
      ang <- if (m > 1) opt_global$par[idx:(idx + len_angles - 1)] else numeric(0)
      if (m > 1) idx <- idx + len_angles
      new_components[[k]]$angles <- ang
      new_components[[k]]$reflection <- ref_vec[k]
      new_components[[k]]$O <- build_O_R(ang, m, ref_vec[k])
    }
    u_pars <- opt_global$par[idx:(idx + len_U_diag + len_U_vec - 1)]
    U_shared <- build_U_R(u_pars, m)
    for (k in seq_len(K)) new_components[[k]]$U <- U_shared

  } else {
    stop("Unknown constraint_mode: ", constraint_mode)
  }

  new_components
}

# Initialization: k-means / hierarchical / random partitioning. Sets up components with empirical mu and U; identity O; moderate gamma and nu.
initialize_mixsit <- function(Y, K, method = "kmeans", constraint_mode = "none") {
  Y <- as.matrix(Y)
  N <- nrow(Y); m <- ncol(Y)

  labels <- switch(
    method,
    # Candidate protocol: exactly one package-native k-means start.
    "kmeans" = stats::kmeans(Y, centers = K, nstart = 5, iter.max = 50)$cluster,
    "hclust" = stats::cutree(stats::hclust(stats::dist(Y), method = "ward.D2"), k = K),
    "random" = sample.int(K, N, replace = TRUE),
    stop("Unknown init method: ", method)
  )

  pi_init <- as.numeric(table(factor(labels, levels = 1:K))) / N
  pi_init[pi_init < 1e-3] <- 1e-3
  pi_init <- pi_init / sum(pi_init)

  components <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k <- which(labels == k)
    if (length(idx_k) < m + 1) idx_k <- sample.int(N, max(m + 1, 5))
    Y_k <- Y[idx_k, , drop = FALSE]

    mu_k <- colMeans(Y_k)
    S_k <- stats::cov(Y_k) + diag(1e-3, m, m)
    U_k <- tryCatch(chol(S_k), error = function(e) diag(sqrt(diag(S_k)), m, m))

    components[[k]] <- list(
      mu = mu_k,
      U = U_k,
      O = diag(m),
      gamma = rep(1.2, m),
      nu = seq(5, 5 + (m - 1) * 0.5, length.out = m),
      angles = if (m > 1) rep(0, m * (m - 1) / 2) else numeric(0),
      reflection = FALSE
    )
  }

  # EV and EE must enter their first E-step inside the constrained model space.
  if (constraint_mode %in% c("U", "full")) {
    centers <- vapply(seq_len(K), function(k) {
      ii <- which(labels == k)
      if (length(ii)) colMeans(Y[ii, , drop = FALSE]) else colMeans(Y)
    }, numeric(m))
    residuals <- Y
    for (i in seq_len(N)) residuals[i, ] <- Y[i, ] - centers[, labels[i]]
    S_pool <- crossprod(residuals) / max(N - K, 1) + diag(1e-3, m, m)
    U_pool <- tryCatch(chol(S_pool),
      error = function(e) diag(sqrt(pmax(diag(S_pool), 1e-8)), m, m))
    for (k in seq_len(K)) components[[k]]$U <- U_pool
  }

  list(pi = pi_init, components = components)
}

# One complete two-cycle AECM iteration. Cycle 1 uses the Z augmentation for pi and gamma.  Cycle 2 recomputes the
# posterior quantities under the intermediate parameters before updating nu and the structural block with the (Z, W) augmentation.  The final observed
# penalized likelihood is checked directly as a numerical safety net.
mixsit_aecm_step <- function(Y, params, expectations_z = NULL,
                             lambda_U = 0, lambda_O = 0,
                             lambda_magnitude = 0, lambda_hetero = 5,
                             constraint_mode = "full",
                             ascent_tol = 1e-10) {
  old_params <- params
  if (is.null(expectations_z))
    expectations_z <- tryCatch(E_Step_cpp(Y, params), error = function(e) NULL)
  if (is.null(expectations_z))
    return(list(params = old_params, expectations = NULL, accepted = FALSE,
                old_llp = -Inf, new_llp = -Inf, reason = "initial E-step failed"))

  old_pen <- compute_penalty_total(old_params, lambda_U, lambda_O,
                                   lambda_magnitude, lambda_hetero)
  old_llp <- expectations_z$loglik - old_pen

  # AECM cycle 1: Z is the missing data and W is integrated out.
  params$pi <- CM_Step_Pi(expectations_z$Z)
  params$components <- CM_Step_Gamma(Y, params, expectations_z, constraint_mode)
  params <- gauge_canonicalise(params, constraint_mode)

  # AECM cycle 2: refresh both memberships and latent-scale moments after the
  # gamma update.  Reusing expectations_z here would not be a valid AECM cycle.
  expectations_w <- tryCatch(E_Step_cpp(Y, params), error = function(e) NULL)
  if (is.null(expectations_w))
    return(list(params = old_params, expectations = expectations_z,
                accepted = FALSE, old_llp = old_llp, new_llp = -Inf,
                reason = "refreshed E-step failed"))
  params$components <- CM_Step_Nu(params, expectations_w,
                                  lambda_magnitude, lambda_hetero,
                                  constraint_mode)
  params$components <- CM_Step_Space(Y, params, expectations_w,
                                     lambda_U, lambda_O, constraint_mode)
  params <- gauge_canonicalise(params, constraint_mode)

  expectations_new <- tryCatch(E_Step_cpp(Y, params), error = function(e) NULL)
  if (is.null(expectations_new))
    return(list(params = old_params, expectations = expectations_z,
                accepted = FALSE, old_llp = old_llp, new_llp = -Inf,
                reason = "final E-step failed"))
  new_pen <- compute_penalty_total(params, lambda_U, lambda_O,
                                   lambda_magnitude, lambda_hetero)
  new_llp <- expectations_new$loglik - new_pen
  slack <- ascent_tol * (1 + abs(old_llp))
  if (!is.finite(new_llp) || new_llp < old_llp - slack) {
    return(list(params = old_params, expectations = expectations_z,
                accepted = FALSE, old_llp = old_llp, new_llp = new_llp,
                reason = "observed penalized likelihood decreased"))
  }
  list(params = params, expectations = expectations_new, accepted = TRUE,
       old_llp = old_llp, new_llp = new_llp, reason = "accepted")
}

# Main fitter: ECM with short-run screening initialization.
#
# Convergence: |ell_p^(t+1) - ell_p^(t)| < tol  OR  iter >= max_iter.
# Returns:
#   params         : list(pi, components)
#   loglik         : trajectory of ell_p across iterations
#   cluster        : hard cluster assignment (length N)
#   posterior      : posterior probability matrix (N x K)
#   converged      : TRUE if stopped at tol; FALSE if the iteration cap is hit or a guarded AECM step is rejected
#   iterations     : number of EM iterations performed
fitMixSIT_ECM <- function(Y, K = 2, constraint_mode = "full",
                          n_short = 10, short_iter = 25,
                          init_methods = c("kmeans", "hclust", "random"),
                          max_iter = 500, tol = 1e-4,
                          tol_type = c("relative", "absolute"),
                          verbose = FALSE,
                          lambda_U = 0, lambda_O = 0,
                          lambda_magnitude = 0, lambda_hetero = 5,
                          flips = TRUE, flip_outer = 6, flip_inner = 60) {

  tol_type <- match.arg(tol_type)

  Y <- as.matrix(Y)
  N <- nrow(Y); m <- ncol(Y)

  # Short-run screening: rank candidates by ell_p
  candidates <- list()
  llp_short <- rep(-Inf, n_short)
  methods_pool <- rep(init_methods, length.out = n_short)

  for (i in seq_len(n_short)) {
    params <- tryCatch(initialize_mixsit(Y, K, methods_pool[i], constraint_mode),
                       error = function(e) NULL)
    if (is.null(params)) next

    expectations <- tryCatch(E_Step_cpp(Y, params), error = function(e) NULL)
    for (iter in seq_len(short_iter)) {
      if (is.null(expectations)) break
      st <- mixsit_aecm_step(
        Y, params, expectations,
        lambda_U, lambda_O, lambda_magnitude, lambda_hetero,
        constraint_mode)
      if (!isTRUE(st$accepted)) break
      params <- st$params
      expectations <- st$expectations
    }
    if (!is.null(expectations)) {
      expectations <- tryCatch(E_Step_cpp(Y, params), error = function(e) NULL)
      if (!is.null(expectations)) {
        pen <- compute_penalty_total(params, lambda_U, lambda_O,
                                     lambda_magnitude, lambda_hetero)
        candidates[[i]] <- params
        llp_short[i] <- expectations$loglik - pen
      }
    }
  }
  best_idx <- which.max(llp_short)
  if (length(best_idx) == 0 || !is.finite(llp_short[best_idx])) {
    stop("Initialization failed: no candidate produced a valid likelihood.")
  }
  params <- candidates[[best_idx]]

  # Main ECM loop with ell_p monitoring
  expectations <- tryCatch(E_Step_cpp(Y, params), error = function(e) NULL)
  if (is.null(expectations)) stop("Initial E-step failed after screening.")
  initial_pen <- compute_penalty_total(params, lambda_U, lambda_O,
                                       lambda_magnitude, lambda_hetero)
  llp_hist <- expectations$loglik - initial_pen
  converged <- FALSE
  iter_used <- 0L
  ascent_rejections <- 0L

  for (iter in seq_len(max_iter)) {
    prev_llp <- tail(llp_hist, 1)
    st <- mixsit_aecm_step(
      Y, params, expectations,
      lambda_U, lambda_O, lambda_magnitude, lambda_hetero,
      constraint_mode)
    iter_used <- iter
    if (!isTRUE(st$accepted)) {
      ascent_rejections <- ascent_rejections + 1L
      if (verbose) cat("AECM step rejected at iter", iter, ":", st$reason, "\n")
      break
    }
    params <- st$params
    expectations <- st$expectations
    curr_llp <- st$new_llp
    llp_hist <- c(llp_hist, curr_llp)
    if (verbose) cat(sprintf("Iter %d: ell_p = %.6f\n", iter, curr_llp))

    delta <- abs(curr_llp - prev_llp)
    criterion <- if (tol_type == "relative")
      delta / (abs(prev_llp) + 1e-10) else delta
    if (is.finite(criterion) && criterion < tol) {
      converged <- TRUE
      break
    }
  }

  # Axis-flip stage  
  # Runs on the converged parameters, alternating guarded flips with ordinary ECM until no flip is accepted. Every candidate is kept only if ell_p rises,
  # so a failure here can only leave `params` where the ECM loop left them.
  n_flips <- 0L
  if (isTRUE(flips)) {
    fp <- tryCatch(mixsit_flip_loop(Y, params, lam_het = lambda_hetero,
                     mode = constraint_mode, outer = flip_outer,
                     inner = flip_inner, lam_U = lambda_U, lam_O = lambda_O,
                     tol = tol), error = function(e) NULL)
    if (!is.null(fp)) { params <- fp$params; n_flips <- fp$n_flips }
  }
  # EE: gamma is free, so identifiability is fixed by one canonical choice per
  # SHARED axis rather than by the parameter space. Likelihood-neutral.
  if (identical(constraint_mode, "full")) {
    cp <- tryCatch(canonicalise_shared_axes(params), error = function(e) NULL)
    if (!is.null(cp)) params <- cp
  }

  # Final E-step and label relabeling by max_j gamma
  expectations <- E_Step_cpp(Y, params)
  if (K > 1) {
    gamma_max <- sapply(params$components, function(comp) max(comp$gamma))
    sort_order <- order(gamma_max)
    params$components <- params$components[sort_order]
    params$pi <- params$pi[sort_order]
    expectations$Z <- expectations$Z[, sort_order]
  }
  final_cluster <- apply(expectations$Z, 1, which.max)

  list(
    params      = params,
    loglik      = llp_hist,
    # ell(Omega) at the final parameters, WITHOUT the shrinkage term. `loglik` above is the ell_p trajectory the algorithm ascends; model-selection
    # criteria must use the raw log-likelihood, so that they mean the same thing for MixSIT as for every competing method.
    loglik_obs  = as.numeric(expectations$loglik),         # final raw ell (NOT ell_p)
    cluster     = final_cluster,
    posterior   = expectations$Z,
    converged   = converged,
    iterations  = iter_used,
    n_flips     = n_flips,
    ascent_rejections = ascent_rejections
  )
}
