# Same directional CM objective and plane rotations as v2.0.0. Evaluate and
# accept steps using nonnegative weighted residual squares, NOT g' A g.
orientation_cm <- function(X,weights,mu,G,max_sweeps=100L,reltol=1e-8) {
  d<-ncol(X);R<-sweep(X,2,mu,'-')
  stopifnot(d>=2L,all(is.finite(R)),all(is.finite(weights)),all(weights>=0),
    all(is.finite(G)),max(abs(crossprod(G)-diag(d)))<1e-6)
  sw<-sqrt(weights)
  project<-function(H)vapply(seq_len(ncol(H)),function(j)as.vector(R%*%H[,j]),numeric(nrow(R)))
  cost<-function(Z)sum((sw*Z)^2)
  Z<-project(G);history<-cost(Z)
  if(!is.finite(history))stop('Nonfinite initial direction objective')
  converged<-FALSE;accepted<-rejected<-0L;termination<-'sweep_cap'
  for(sweep in seq_len(max_sweeps)) {
    G_before<-G;Z_before<-Z;before<-tail(history,1);sweep_accepted<-0L
    for(j in seq_len(d-1L))for(k in (j+1L):d) {
      uj<-sw[,j]*Z[,j];vj<-sw[,j]*Z[,k]
      uk<-sw[,k]*Z[,j];vk<-sw[,k]*Z[,k]
      scale<-max(abs(c(uj,vj,uk,vk)))
      if(scale==0)next
      if(!is.finite(scale)){rejected<-rejected+1L;next}
      # Scale the coefficients BEFORE squaring; their angle is scale invariant.
      uj<-uj/scale;vj<-vj/scale;uk<-uk/scale;vk<-vk/scale
      a<-sum(uj^2+vk^2);b<-sum(vj^2+uk^2);c<-sum(uj*vj-uk*vk)
      h<-(a-b)/2;radius<-sqrt(h*h+c*c)
      if(radius==0)next
      # Stable form of radius+h when h<0 (small improvement / large b).
      improvement<-if(h<0)c*c/(radius-h) else radius+h
      if(!is.finite(improvement)||improvement<=0)next
      angle<-.5*atan2(-c,-h);cs<-cos(angle);sn<-sin(angle)
      u<-G[,j];v<-G[,k];proposal<-cbind(cs*u+sn*v,-sn*u+cs*v)
      Znew<-project(proposal)
      oldcost<-sum((sw[,j]*Z[,j])^2+(sw[,k]*Z[,k])^2)
      newcost<-sum((sw[,j]*Znew[,1])^2+(sw[,k]*Znew[,2])^2)
      # Strict acceptance: no tolerated increase. No threshold relaxation.
      if(is.finite(newcost)&&newcost<oldcost) {
        G[,j]<-proposal[,1];G[,k]<-proposal[,2]
        Z[,j]<-Znew[,1];Z[,k]<-Znew[,2]
        accepted<-accepted+1L;sweep_accepted<-sweep_accepted+1L
      } else rejected<-rejected+1L
    }
    value<-cost(Z)
    if(!is.finite(value))stop('Nonfinite direction objective')
    if(value>before) {
      # Full objective is a second independent acceptance check. Restore the
      # last feasible state if summation roundoff defeats the pair checks.
      G<-G_before;Z<-Z_before;value<-before;termination<-'numerical_stagnation'
      history<-c(history,value);break
    }
    history<-c(history,value)
    if(before-value<=reltol*(1+abs(before))) {
      converged<-sweep_accepted>0L;termination<-if(converged)'relative_objective' else 'no_accepted_rotation'
      break
    }
  }
  if(max(abs(crossprod(G)-diag(d)))>1e-7)stop('Direction CM lost orthogonality')
  list(G=G,objective=history,sweeps=sweep,converged=converged,
       termination=termination,accepted_rotations=accepted,rejected_rotations=rejected)
}
gamma_stable <- function(X,z=NULL,w=NULL,mu,Lambda,Gamma0,method='BFGS',
                         iter.max=5000,reltol=1e-15,trace=0) {
  if(is.null(z))z<-rep(1,nrow(X));if(is.null(w))w<-matrix(1,nrow(X),ncol(X))
  orientation_cm(X,sweep(w*z,2,diag(Lambda),'/'),mu,Gamma0,
    max_sweeps=iter.max,reltol=reltol)$G
}
