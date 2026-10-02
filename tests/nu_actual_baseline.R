# Regression: CM_Step_Nu must compare against the actual incoming nu,
# even when the optimizer's initialization floors change its encoded start.
args0 <- commandArgs(FALSE)
this <- sub('^--file=', '', args0[grep('^--file=', args0)])
root <- Sys.getenv('MIXSIT_TEST_ROOT', unset = '')
if (!nzchar(root)) root <- normalizePath(file.path(dirname(this), '..'), winslash = '/')
.libPaths(c(file.path(root, '.r-lib'), .libPaths()))
source(file.path(root, 'Stage1', 'MixSIT_rcpp.R'))
objective <- function(nu, ex, lm, lh) {
  h <- nu / 2; lv <- log(nu); z <- ex$Z[, 1]
  td <- vapply(seq_along(nu), function(j)
    sum(z * (ex$Kappa[, j, 1] - ex$W[, j, 1])), numeric(1))
  -(sum(z) * sum(h * log(h) - lgamma(h)) + sum(h * td)) +
    lm * sum((lv - log(10))^2) + lh * sum((lv - mean(lv))^2)
}
make_params <- function(nu) list(pi = 1, components = list(list(
  mu = c(0, 0), U = diag(2), O = diag(2), gamma = c(1.2, 1.2),
  nu = nu, angles = 0, reflection = FALSE)))
for (lh in c(0, 5)) {
  nu <- c(5, 5.01); Y <- matrix(0, 100, 2)
  for (j in 1:2) {
    base <- qt((seq_len(100) - .5) / 100, df = nu[j])
    target <- 2 * lh * (log(nu[j]) - mean(log(nu))) / nu[j]
    score <- function(s) {
      x <- base * s; v <- nu[j]
      sum(.5 * (digamma((v + 1)/2) - digamma(v/2) - 1/v) -
          .5 * log1p(x^2/v) + .5 * (v+1) * x^2/(v*(v+x^2))) - target
    }
    x <- base * uniroot(score, c(.01, 100), tol = 1e-12)$root
    Y[, j] <- x * 1.2^sign(x)
  }
  p <- make_params(nu); ex <- E_Step_cpp(Y, p)
  a <- CM_Step_Nu(p, ex, 0, lh, 'U')
  stopifnot(identical(a[[1]]$nu, nu))
  for (mode in c('U', 'full'))
    for (v in list(c(1.2, 1.21), c(5, 5), c(5, 5+1e-10), c(10, 15))) {
      p <- make_params(v); ex <- E_Step_cpp(Y, p)
      for (lm in c(0, 2)) {
        a <- CM_Step_Nu(p, ex, lm, lh, mode)
        stopifnot(objective(a[[1]]$nu, ex, lm, lh) <= objective(v, ex, lm, lh))
        if (identical(v, c(10, 15)) && lm == 0)
          stopifnot(objective(a[[1]]$nu, ex, lm, lh) < objective(v, ex, lm, lh))
      }
    }
}
original_solver <- .cm2_solve
.cm2_solve <- function(par0, vg, fd_fn, control)
  list(par = rep(Inf, length(par0)), objective = Inf)
for (mode in c('U', 'full')) {
  a <- CM_Step_Nu(p, ex, 0, 5, mode)
  stopifnot(identical(a[[1]]$nu, p$components[[1]]$nu))
}
.cm2_solve <- original_solver
cat('actual incoming nu baseline regression tests passed\n')
