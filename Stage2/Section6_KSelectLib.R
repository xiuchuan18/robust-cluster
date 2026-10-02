# Section6_KSelectLib.R
#
# Truths and information criteria for the Section 6.2 study of choosing the number of components.

SEP_K <- 5      # default class separation, in the units the generators use

# `sep` is an argument rather than a constant so that the same truths can be run
# down a difficulty ladder. The main grid uses SEP_K = 5, where the components
# are unambiguous and a wrong K can only mean the model needed a different NUMBER
# of components for the shape. Lowering it puts the ordinary difficulty of
# finding clusters back in, which is where any selection rule eventually breaks
# and is worth locating rather than avoiding.

gen_kcell <- function(cell, N, seed, sep = SEP_K) {
  set.seed(seed)
  switch(cell,

    G3 = { d <- 5; K <- 3
      mu <- list(rep(0, d), mu_at(d, 1:2, sep), mu_at(d, 3:4, sep))
      nk <- as.vector(rmultinom(1, N, rep(1/K, K)))
      Y <- do.call(rbind, lapply(1:K, function(k)
        mvtnorm::rmvnorm(nk[k], mu[[k]], diag(d))))
      list(Y = Y, labels = rep(1:K, nk), K = K) },

    B1 = { d <- 5; K <- 1
      O <- qr.Q(qr(matrix(rnorm(d^2), d, d)))
      Y <- assemble(N, rep(0, d), diag(d), O,
                    seq(4, 9, length.out = d), seq(0.55, 2.1, length.out = d))
      list(Y = Y, labels = rep(1L, N), K = K) },

    # The K = 1 pair. B1 is drawn from the MixSIT family itself, so MixSIT is correctly specified there and "it returns K = 1" is close to a tautology;
    # it is kept only as the aligned reference. U1 is the one that carries the claim: a single unrestricted skew-t population, which NO fitted model can
    # represent exactly, so returning K = 1 there is a real result rather than a home fixture. Neither cell contributes to any ARI table -- with one true
    # class the index is not defined -- and both are reported as selection results only.
    U1 = { d <- 5; K <- 1
      Y <- rmst_sdb(N, rep(0, d), get_fixed_covariance(d, 1),
                    c(4, 3, 2, -3, 2), 20)
      list(Y = Y, labels = rep(1L, N), K = K) },

    B2 = { d <- 5; K <- 2
      mu <- list(rep(0, d), mu_at(d, 1:2, sep))
      gm <- list(seq(0.55, 2.1, length.out = d), seq(1.8, 0.6, length.out = d))
      nu <- list(seq(4, 9, length.out = d), seq(3, 7, length.out = d))
      O  <- lapply(1:K, function(k) qr.Q(qr(matrix(rnorm(d^2), d, d))))
      nk <- as.vector(rmultinom(1, N, rep(1/K, K)))
      Y <- do.call(rbind, lapply(1:K, function(k)
        assemble(nk[k], mu[[k]], diag(d), O[[k]], nu[[k]], gm[[k]])))
      list(Y = Y, labels = rep(1:K, nk), K = K) },

    A3 = { d <- 10; K <- 3
      mu <- list(rep(0, d), mu_at(d, 1:3, sep), mu_at(d, 4:6, sep))
      gm <- list(seq(0.6, 1.8, length.out = d), seq(1.25, 2.25, length.out = d),
                 seq(0.5, 1.5, length.out = d))
      nu <- list(seq(4, 8, length.out = d), seq(3, 6, length.out = d),
                 seq(5, 10, length.out = d))
      O  <- lapply(1:K, function(k) qr.Q(qr(matrix(rnorm(d^2), d, d))))
      nk <- as.vector(rmultinom(1, N, rep(1/K, K)))
      Y <- do.call(rbind, lapply(1:K, function(k)
        assemble(nk[k], mu[[k]], diag(d), O[[k]], nu[[k]], gm[[k]])))
      list(Y = Y, labels = rep(1:K, nk), K = K) },

    U3 = { d <- 5; K <- 3
      mu <- list(rep(0, d), mu_at(d, 1:2, sep), mu_at(d, 3:4, sep))
      dl <- list(c(4, 3, 2, -3, 2), c(5, 2, 3, -2, 3), c(4, -2, -2, -2, 4))
      sg <- lapply(1:3, function(k) get_fixed_covariance(d, k))
      nk <- as.vector(rmultinom(1, N, rep(1/K, K)))
      Y <- do.call(rbind, lapply(1:K, function(k)
        rmst_sdb(nk[k], mu[[k]], sg[[k]], dl[[k]], 20)))
      list(Y = Y, labels = rep(1:K, nk), K = K) },

    T4 = { d <- 5; K <- 4
      mu <- list(rep(0, d), mu_at(d, 1:2, sep), mu_at(d, 3:4, sep),
                 mu_at(d, c(1, 3), -sep))
      nk <- as.vector(rmultinom(1, N, rep(1/K, K)))
      Y <- do.call(rbind, lapply(1:K, function(k) {
        G <- rgamma(nk[k], 2, 2)                       # t_4 scale mixture
        sweep(mvtnorm::rmvnorm(nk[k], rep(0, d), diag(d)) / sqrt(G), 2,
              mu[[k]], "+") }))
      list(Y = Y, labels = rep(1:K, nk), K = K) },

    C3 = { d <- 10; K <- 3
      mu <- list(rep(0, d), mu_at(d, 1:4, sep / 2), mu_at(d, 5:8, sep / 2))
      nk <- as.vector(rmultinom(1, N, rep(1/K, K)))
      Y <- do.call(rbind, lapply(1:K, function(k)
        mvtnorm::rmvnorm(nk[k], mu[[k]], diag(d))))
      C <- matrix(runif(N * d) < 0.05, N, d)
      Y[C] <- Y[C] + rnorm(sum(C), 0, 10)
      list(Y = Y, labels = rep(1:K, nk), K = K) },

    stop("unknown K-selection cell: ", cell))
}

KCELLS <- c("G3", "B1", "U1", "B2", "A3", "U3", "T4", "C3")

# criteria, uniformly
npar_of <- function(method, K, d) count_parameters(method, K, d)

entropy_of <- function(z, n) posterior_entropy(z, n)

ic_of <- function(loglik, z, method, K, n, d) {
  p <- npar_of(method, K, d)
  bic <- -2 * loglik + p * log(n)
  ent <- entropy_of(z, n)
  list(npar = p, AIC = -2 * loglik + 2 * p, BIC = bic,
       ICL = if (is.na(ent)) NA_real_ else bic + 2 * ent)
}
