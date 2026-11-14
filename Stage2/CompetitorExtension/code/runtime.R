repo<-normalizePath('../..',mustWork=TRUE)
.libPaths(c(file.path(repo,'.r-lib'),.libPaths()))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
source('code/00_RNG.R');rng_fix()
PORT_VERSION <- 'v2.0.1-direct-residual-jacobi'
logadd <- function(a,b) { m <- pmax(a,b); m+log(exp(a-m)+exp(b-m)) }
row_lse <- function(x) { m <- apply(x,1,max); m+log(rowSums(exp(x-m))) }
atomic_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(tmpdir=dirname(path)); saveRDS(x,tmp)
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
# Algebraically identical to MSclust::Gamma.update's trace objective; no
# change to coordinates, numerical derivatives, optim method or tolerances.
gamma_fast <- function(X,z=NULL,w=NULL,mu,Lambda,Gamma0,method="BFGS",
                       iter.max=5000,reltol=1e-15,trace=0) {
  n <- nrow(X); d <- ncol(X)
  if(is.null(z)) z <- rep(1,n)
  if(is.null(w)) w <- matrix(1,n,d)
  fac <- MSclust:::QPl(Gamma0)
  xc <- sweep(X,2,mu,"-"); weights <- sweep(w*z,2,diag(Lambda),"/")
  f <- function(par) sum(weights*(xc%*%MSclust:::PlQ(fac$P,par))^2)
  op <- optim(fac$l,f,method=method,control=list(maxit=iter.max,reltol=reltol,trace=trace))
  MSclust:::PlQ(fac$P,op$par)
}


