source('code/runtime.R');source('code/orientation.R');source('code/tortora2025.R')
fit_new <- function(code,Y,K,seed) {
  rng_fix(); set.seed(seed)
  if(code=="FW2014_MSt") {
    km <- stats::kmeans(Y,K,nstart=5)
    z <- mclust::unmap(km$cluster)
    f <- MSclust::mst
    e <- new.env(parent=environment(f)); e$Gamma.update <- gamma_fast
    environment(f) <- e
    fit <- f(Y,K,ini="manual",sz=z)
    list(fit=fit,z=fit$z,cluster=fit$cluster,loglik=fit$loglik,
         iter=fit$iter.stop,converged=NA,at_cap=fit$iter.stop>=200,
         termination=if(fit$iter.stop>=200)"iteration_cap" else "native_stop_not_certified",
         parameter_count=fit$npar,km=km$cluster)
  } else if(code=="ZH2019_ICA") {
    # Use the published native compiled engine; replace only the R wrapper's
    # fixed seed and single-start kmeans. Avoid its fixed-name tempfile race.
    km <- stats::kmeans(Y,K,nstart=5)
    z <- matrix(0,nrow(Y),K); z[cbind(seq_len(nrow(Y)),km$cluster)] <- 1
    loadNamespace("icamix")
    fit <- .Call("icamix_EMInterwovenFastICA",PACKAGE="icamix",t(Y),z,
                 0,301L,150L,1e-6,FALSE,TRUE)
    it <- nrow(fit$Lambdas)-1L
    delta <- max(abs(fit$Lambdas[it+1,]-fit$Lambdas[it,]))
    list(fit=fit,z=fit$MembershipProbs,cluster=max.col(fit$MembershipProbs,ties.method="first"),
         loglik=as.numeric(tail(fit$ObjValue,1)),iter=it,converged=delta<=1e-6,
         at_cap=it>=300,termination=if(delta<=1e-6)"native_weight_tolerance" else "iteration_cap",
         parameter_count=NA_real_,km=km$cluster)
  } else if(code %in% c("T2025_Manly","T2025_Power")) {
    km <- stats::kmeans(Y,K,nstart=5)
    fit <- fit_t25(Y,K,if(code=="T2025_Manly")"manly" else "power",km)
    list(fit=fit,z=fit$z,cluster=fit$cluster,loglik=fit$loglik,iter=fit$iter,
         converged=fit$converged,at_cap=fit$iter>=500,termination=fit$termination,
         parameter_count=K-1+K*(4*ncol(Y)+ncol(Y)*(ncol(Y)+1)/2),km=km$cluster)
  } else stop("Unsupported model: ",code)
}
