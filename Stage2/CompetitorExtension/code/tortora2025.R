# Independent ECM implementation of Tortora, Punzo & Franczak (2025),
# doi:10.1007/s00362-025-01723-9, equations 7,9,16--22 and CM step 4.
# Initialization is the USER'S existing kmeans(nstart=5) protocol, not paper PAM.
# Unspecified numerical defaults: 500 ECM iterations, Aitken tolerance 1e-5;
# NM transform optimization uses R optim defaults (maxit=500).
# v2: directional CM uses at most 500 sweeps of exact plane rotations;
# the objective is unchanged, with relative decrease tolerance 1e-5.
transform_t25 <- function(X,lambda,type) {
  n <- nrow(X); d <- ncol(X); T <- X; J <- matrix(0,n,d)
  for(j in seq_len(d)) {
    x <- X[,j]; a <- lambda[j]
    if(type=="manly") {
      T[,j] <- if(abs(a)<1e-10) x else expm1(a*x)/a
      J[,j] <- a*x
    } else {
      pos <- x>=0; lp <- log1p(x[pos]); ln <- log1p(-x[!pos])
      T[pos,j] <- if(abs(a)<1e-10) lp else expm1(a*lp)/a
      b <- 2-a
      T[!pos,j] <- if(abs(b)<1e-10) -ln else -expm1(b*ln)/b
      J[,j] <- (a-1)*sign(x)*log1p(abs(x))
    }
  }
  list(T=T,logJ=rowSums(J))
}
t25_estep <- function(X,pars,type) {
  n <- nrow(X); K <- length(pars); d <- ncol(X)
  ld <- matrix(0,n,K); v <- vector("list",K)
  for(k in seq_len(K)) {
    p <- pars[[k]]; tr <- transform_t25(X,p$lambda,type)
    r <- sweep(tr$T,2,p$mu,"-")%*%p$G
    good <- sweep(-.5*sweep(r^2,2,p$L,"/"),2,
                  log(p$alpha)-.5*log(2*pi*p$L),"+")
    bad <- sweep(-.5*sweep(r^2,2,p$L*p$eta,"/"),2,
                 log1p(-p$alpha)-.5*log(2*pi*p$L*p$eta),"+")
    den <- logadd(good,bad); v[[k]] <- exp(good-den)
    ld[,k] <- log(p$pi)+rowSums(den)+tr$logJ
  }
  lse <- row_lse(ld)
  list(z=exp(ld-lse),v=v,loglik=sum(lse))
}
fit_t25 <- function(X,K,type,km,maxiter=500L,tol=1e-5) {
  n <- nrow(X); d <- ncol(X)
  pars <- lapply(seq_len(K),function(k) {
    xx <- X[km$cluster==k,,drop=FALSE]; mu <- colMeans(xx)
    S <- crossprod(sweep(xx,2,mu,"-"))/nrow(xx); ee <- eigen(S,symmetric=TRUE)
    if(any(ee$values<=0)) stop("Singular initial within-cluster covariance")
    list(pi=nrow(xx)/n,mu=mu,G=ee$vectors,L=ee$values,
         alpha=rep(.99,d),eta=rep(1.01,d),lambda=rep(if(type=="manly")0 else 1,d))
  })
  es <- t25_estep(X,pars,type); ll <- es$loglik; converged <- FALSE
  termination <- "iteration_cap"; optcodes <- list()
  for(it in seq_len(maxiter)) {
    old <- pars
    for(k in seq_len(K)) {
      p <- pars[[k]]; z <- es$z[,k]; v <- es$v[[k]]; nk <- sum(z)
      if(!is.finite(nk)||nk<1e-10) stop("Empty component")
      tr <- transform_t25(X,p$lambda,type); T <- tr$T
      p$pi <- nk/n
      p$alpha <- pmin(1-1e-10,pmax(.5,colSums(v*z)/nk))
      w <- v+sweep(1-v,2,p$eta,"/")
      zw <- w*z
      # Eq 20: weighted means in the current principal-coordinate space.
      p$mu <- as.vector(p$G%*%(colSums(zw*(T%*%p$G))/colSums(zw)))
      r <- sweep(T,2,p$mu,"-")%*%p$G
      p$L <- colSums(zw*r^2)/nk
      if(any(!is.finite(p$L))||any(p$L<=0)) stop("Degenerate component scale")
      badw <- (1-v)*z
      p$eta <- pmax(1.001,colSums(badw*r^2)/colSums(badw)/p$L)
      if(any(!is.finite(p$eta))) stop("Undefined contamination inflation")
      w <- v+sweep(1-v,2,p$eta,"/")
      p$G <- gamma_stable(T,z,w,p$mu,diag(p$L),p$G,iter.max=500,reltol=1e-5)
      weights <- sweep(w*z,2,p$L,"/")
      fn <- function(lambda) {
        a <- transform_t25(X,lambda,type)
        if(any(!is.finite(a$T))||any(!is.finite(a$logJ))) return(.Machine$double.xmax/100)
        r <- sweep(a$T,2,p$mu,"-")%*%p$G
        val <- .5*sum(weights*r^2)-sum(z*a$logJ)
        if(is.finite(val)) val else .Machine$double.xmax/100
      }
      op <- optim(p$lambda,fn,method="Nelder-Mead")
      optcodes[[length(optcodes)+1]] <- op$convergence
      if(op$value<=fn(p$lambda)) p$lambda <- op$par
      pars[[k]] <- p
    }
    nextes <- t25_estep(X,pars,type)
    if(!is.finite(nextes$loglik)) stop("Nonfinite observed likelihood")
    if(nextes$loglik < tail(ll,1)-1e-7*(1+abs(tail(ll,1)))) {
      pars <- old; termination <- "likelihood_decrease"; break
    }
    es <- nextes; ll <- c(ll,es$loglik)
    if(length(ll)>=3) {
      vll <- tail(ll,3); denom <- vll[2]-vll[1]
      a <- (vll[3]-vll[2])/denom
      gap <- (vll[3]-vll[2])/(1-a)
      if(is.finite(a)&&a>=0&&a<1&&is.finite(gap)&&gap>=0&&gap<tol) {
        converged <- TRUE; termination <- "aitken"; break
      }
    }
    if(it%%10==0) {cat("iteration",it,"loglik",tail(ll,1),"\n"); flush.console()}
  }
  list(params=pars,z=es$z,cluster=max.col(es$z,ties.method="first"),
       loglik=tail(ll,1),trajectory=ll,iter=it,converged=converged,
       termination=termination,optim_codes=unlist(optcodes),
       implementation="independent_ECM_paper_equations_Jacobi_direction_v2")
}
