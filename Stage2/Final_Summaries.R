# Numerical publication summaries from aggregate results; run from Stage2/.
read_final <- function(stem) readRDS(file.path("results", paste0(stem, ".rds")))
save_summary <- function(x, stem) {
  write.csv(x, file.path("results", paste0(stem, ".csv")), row.names=FALSE)
}
rmse <- function(x) sqrt(mean(x^2, na.rm=TRUE))
perfit <- function(d, col) unique(d[, c("rep", col)])[[col]]
percomp <- function(d, col) unique(d[, c("rep", "k", col)])[[col]]
angle <- function(d) mean(acos(pmin(d$cos_abs, 1)) * 180 / pi)

R <- read_final("Section6_Recovery")
RS <- do.call(rbind, lapply(split(R, list(R$variant, R$N), drop=TRUE), function(d)
  data.frame(variant=d$variant[1], N=d$N[1],
    ARI=mean(perfit(d,"ari")), pi=mean(perfit(d,"pi_err")),
    U=mean(perfit(d,"U_err")), mu=rmse(percomp(d,"mu_err")),
    lg=rmse(d$lg_err), lnu=rmse(d$lnu_err), angle=angle(d))))
RS <- RS[order(RS$variant, RS$N), ]
save_summary(RS, "Section6_Recovery_Summary")
slopes <- do.call(rbind, lapply(split(RS, RS$variant), function(d) {
  hi <- d$N >= median(d$N)
  v <- vapply(c("pi", "U", "mu", "lg", "lnu", "angle"), function(nm) {
    ok <- hi & is.finite(d[[nm]]) & d[[nm]] > 0
    if (sum(ok) < 2) return(NA_real_)
    unname(coef(lm(log(d[[nm]][ok]) ~ log(d$N[ok])))[2])
  }, 0)
  data.frame(variant=d$variant[1], as.list(v), check.names=FALSE)
}))
save_summary(slopes, "Section6_Recovery_Slopes")

L <- read_final("Section6_LambdaSensitivity")
LF <- do.call(rbind, lapply(split(L, list(L$lam,L$N,L$rep), drop=TRUE), function(d)
  data.frame(lam=d$lam[1], N=d$N[1], rep=d$rep[1], ARI=d$ari[1],
    loglik=d$loglik[1], bic=d$bic[1], lg=rmse(d$lg_err), lnu=rmse(d$lnu_err),
    mu=rmse(d$mu_err[!duplicated(d$k)]), U=d$U_err[1], angle=angle(d))))
LS <- do.call(rbind, lapply(split(LF, LF$N), function(d) {
  lams <- sort(unique(d$lam))
  paired <- Reduce(intersect, lapply(lams, function(l) d$rep[d$lam==l]))
  d <- d[d$rep %in% paired, ]
  do.call(rbind, lapply(lams, function(l) {
    z <- d[d$lam==l, ]
    data.frame(N=d$N[1], lam=l, paired_reps=length(paired),
      as.list(vapply(c("ARI","loglik","bic","lg","lnu","mu","U","angle"),
                     function(nm) mean(z[[nm]]), 0)), check.names=FALSE)
  }))
}))
save_summary(LS, "Section6_Lambda_Summary")

KS <- read_final("Section6_KSelect")
grid <- sort(unique(KS$K))
groups <- split(KS, list(KS$cell,KS$method,KS$rep), drop=TRUE)
selection <- list()
for (crit in c("AIC","BIC","ICL")) {
  chosen <- do.call(rbind, lapply(groups, function(z) {
    z <- z[order(z$K), ]
    if (nrow(z)!=length(grid) || !identical(z$K,grid) ||
        !all(z$ok %in% TRUE) || !all(is.finite(z[[crit]]))) return(NULL)
    data.frame(cell=z$cell[1], method=z$method[1], rep=z$rep[1],
               K_true=z$K_true[1], K_hat=z$K[which.min(z[[crit]])])
  }))
  cells <- do.call(rbind, lapply(split(KS,list(KS$cell,KS$method),drop=TRUE), function(z) {
    g <- chosen[chosen$cell==z$cell[1] & chosen$method==z$method[1], ]
    data.frame(scope="cell", criterion=crit, cell=z$cell[1], method=z$method[1],
      cells=1L, complete_paths=nrow(g), attempted_paths=length(unique(z$rep)),
      K_true=z$K_true[1], mean_K=if(nrow(g)) mean(g$K_hat) else NA_real_,
      hit=if(nrow(g)) mean(g$K_hat==g$K_true) else NA_real_,
      error=if(nrow(g)) mean(g$K_hat-g$K_true) else NA_real_)
  }))
  overall <- do.call(rbind, lapply(split(cells,cells$method), function(d)
    data.frame(scope="equal_scenario", criterion=crit, cell="ALL", method=d$method[1],
      cells=sum(d$complete_paths>0), complete_paths=sum(d$complete_paths),
      attempted_paths=sum(d$attempted_paths), K_true=NA_real_,
      mean_K=mean(d$mean_K,na.rm=TRUE), hit=mean(d$hit,na.rm=TRUE),
      error=mean(d$error,na.rm=TRUE))))
  selection[[crit]] <- rbind(cells,overall)
}
save_summary(do.call(rbind,selection), "Section6_KSelect_Summary")

S <- read_final("Section6_Scale")
S <- S[S$ok %in% c(TRUE,NA) & is.finite(S$t_tot), ]
sw <- list(N=S[S$m==10 & S$K==3, ], m=S[S$N==800 & S$K==3, ],
           K=S[S$N==800 & S$m==10, ])
expo <- function(d,x,y) {
  a <- aggregate(d[[y]], list(x=d[[x]]), function(z) mean(z,na.rm=TRUE))
  a <- a[is.finite(a$x) & a$x>0 & is.finite(a[[2]]) & a[[2]]>0, ]
  if(nrow(a)<3) return(NA_real_)
  unname(coef(lm(log(a[[2]]) ~ log(a$x)))[2])
}
SS <- do.call(rbind,lapply(c("GMM","MSGHD","MixSIT_None","MixSIT_U","MixSIT_Full"),function(m) {
  d <- S[S$method==m, ]; dm <- sw$m[sw$m$method==m, ]
  dm$periter <- dm$t_base / dm$iters
  ex <- vapply(names(sw), function(x) expo(sw[[x]][sw[[x]]$method==m, ],x,"t_tot"),0)
  data.frame(method=m, N_exponent=ex[["N"]], m_exponent=ex[["m"]], K_exponent=ex[["K"]],
    periter_m=if(all(is.na(dm$iters))) NA_real_ else expo(dm,"m","periter"),
    iterations=if(all(is.na(d$iters))) NA_real_ else mean(d$iters,na.rm=TRUE),
    atcap=if(all(is.na(d$iters))) NA_real_ else mean(d$iters>=if(m=="MSGHD")100 else 500,na.rm=TRUE),
    anchor=mean(d$t_tot[d$N==800 & d$m==10 & d$K==3]),
    m20=mean(d$t_tot[d$N==800 & d$m==20 & d$K==3]))
}))
save_summary(SS, "Section6_Scale_Summary")
cat("wrote five formal numerical summaries from aggregate results\n")
