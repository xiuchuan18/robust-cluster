# Log-nu heterogeneity penalty for the formal lambda-sensitivity study.
# Install after sourcing Stage1. lambda=0 disables this penalty.
PEN <- new.env(parent = emptyenv())
PEN$type <- "log"; PEN$lambda <- 5
install_penalty <- function(type = "log", lambda = 5) {
  stopifnot(identical(type, "log"), length(lambda) == 1L,
            is.finite(lambda), lambda >= 0)
  PEN$type <- type; PEN$lambda <- lambda
  pen_val <- function(nu) {
    if (PEN$lambda <= 0) return(0)
    l <- log(nu); PEN$lambda * sum((l - mean(l))^2)
  }
  pen_grad <- function(nu) {
    if (PEN$lambda <= 0) return(numeric(length(nu)))
    l <- log(nu); 2 * PEN$lambda * (l - mean(l)) / nu
  }
  PEN$val <- pen_val; PEN$grad <- pen_grad
  if (exists("cm2_coord_cpp", envir = globalenv(), inherits = FALSE))
    rm("cm2_coord_cpp", envir = globalenv())
  assign("compute_penalty_total", function(params, lambda_U, lambda_O,
                                           lambda_magnitude, lambda_hetero) {
    sum(vapply(params$components, function(cmp) PEN$val(cmp$nu) + 0, 0))
  }, envir = globalenv())
  assign("CM_Step_Nu", function(params, expectations,
                                lambda_magnitude, lambda_hetero,
                                constraint_mode = "U") {
    K <- length(params$pi); m <- length(params$components[[1]]$nu)
    Z_hat <- expectations$Z; W_hat <- expectations$W
    Kappa_hat <- expectations$Kappa
    new_components <- params$components

    for (k in seq_len(K)) {
      theta <- new_components[[k]]
      z_k <- Z_hat[, k]; sum_z <- sum(z_k)
      if (sum_z < 1e-6) next

      tilde_d <- vapply(seq_len(m), function(j)
        sum(z_k * (Kappa_hat[, j, k] - W_hat[, j, k])), 0)

      if (identical(constraint_mode, "full")) {
        eta0 <- log(pmax(theta$nu, 1e-8))
        vg <- local({
          cp <- NULL; cr <- NULL
          ev <- function(eta) {
            if (is.null(cp) || !identical(eta, cp)) {
              nu <- exp(eta); v <- nu / 2
              val <- -(sum_z * sum(v * log(v) - lgamma(v)) +
                       sum(v * tilde_d)) + PEN$val(nu)
              gnu <- -(sum_z / 2) * (log(v) + 1 - digamma(v)) -
                      tilde_d / 2 + PEN$grad(nu)
              cr <<- list(value = val, grad = gnu * nu); cp <<- eta
            }
            cr
          }
          list(fn = function(eta) ev(eta)$value,
               gr = function(eta) ev(eta)$grad)
        })
        fd <- function(eta) {
          nu <- exp(eta); v <- nu / 2
          -(sum_z * sum(v * log(v) - lgamma(v)) +
            sum(v * tilde_d)) + PEN$val(nu)
        }
        opt_nu <- .cm2_solve(eta0, vg, fd,
                             list(iter.max = 100, eval.max = 100,
                                  rel.tol = 1e-6))
        new_components[[k]]$nu <- exp(opt_nu$par)
      } else {
        vg <- local({
          cp <- NULL; cr <- NULL
          ev <- function(dl) {
            if (is.null(cp) || !identical(dl, cp)) {
              nu <- cumsum(exp(dl)); v <- nu / 2
              val <- -(sum_z * sum(v * log(v) - lgamma(v)) +
                       sum(v * tilde_d)) + PEN$val(nu)
              gnu <- -(sum_z / 2) * (log(v) + 1 - digamma(v)) -
                      tilde_d / 2 + PEN$grad(nu)
              cr <<- list(value = val,
                           grad = exp(dl) * rev(cumsum(rev(gnu))))
              cp <<- dl
            }
            cr
          }
          list(fn = function(dl) ev(dl)$value,
               gr = function(dl) ev(dl)$grad)
        })
        fd <- function(dl) {
          nu <- cumsum(exp(dl)); v <- nu / 2
          -(sum_z * sum(v * log(v) - lgamma(v)) +
            sum(v * tilde_d)) + PEN$val(nu)
        }
        sorted_nu <- sort(theta$nu)
        delta_init <- numeric(m)
        delta_init[1] <- log(max(sorted_nu[1], 2))
        if (m > 1) delta_init[2:m] <- log(pmax(diff(sorted_nu), 0.1))
        opt_nu <- .cm2_solve(delta_init, vg, fd,
                             list(iter.max = 100, eval.max = 100,
                                  rel.tol = 1e-6))
        new_components[[k]]$nu <- cumsum(exp(opt_nu$par))
      }
      # The installed penalty must guard the actual incoming nu, not a potentially floored encoded start, just like the production update.
      objective_actual_nu <- function(nu) {
        if (any(!is.finite(nu)) || any(nu <= 0)) return(Inf)
        v <- nu / 2
        -(sum_z * sum(v * log(v) - lgamma(v)) +
          sum(v * tilde_d)) + PEN$val(nu)
      }
      old_objective <- objective_actual_nu(theta$nu)
      new_objective <- objective_actual_nu(new_components[[k]]$nu)
      if (!is.finite(old_objective) || !is.finite(new_objective) ||
          new_objective > old_objective)
        new_components[[k]]$nu <- theta$nu
    }
    new_components
  }, envir = globalenv())
  invisible(TRUE)
}

check_penalty_gradient <- function(m = 10, tol = 1e-6) {
  set.seed(1)
  nu <- sort(runif(m, 2.5, 12))
  num <- vapply(seq_len(m), function(j) {
    h <- 1e-6; a <- nu; b <- nu; a[j] <- a[j] + h; b[j] <- b[j] - h
    (PEN$val(a) - PEN$val(b)) / (2 * h) }, 0)
  max(abs(num - PEN$grad(nu))) < tol * max(1, max(abs(num)))
}
