repo<-normalizePath('../..');stage<-file.path(repo,'Stage2');.libPaths(c(file.path(repo,'.r-lib'),.libPaths()));base<-'replay/section6'
jobs<-read.csv('../results/competitor-extension/section6.csv')
library(clue)
z<-do.call(rbind,lapply(seq_len(nrow(jobs)),function(i){
 j<-jobs[i,];o<-readRDS(file.path(base,'fits',sprintf('%s_K%d.rds',j$id,j$K)));q<-o$row
 if(j$section=='recovery'){
 x<-readRDS(file.path(base,'data',paste0(j$id,'.rds')));a<-o$fit
 perm<-as.integer(solve_LSAP(table(factor(x$labels,levels=1:3),factor(a$cluster,levels=1:3)),maximum=TRUE))
 normcol<-function(m)sweep(m,2,sqrt(colSums(m*m)),'/')
 angles<-sapply(1:3,function(g){
 k<-perm[g];w<-matrix(a$fit$WMtrs[[k]],j$d);sig<-matrix(a$fit$OriginalSignals[[k]],j$d)
 dif<-w%*%t(x$Y)-sig; stopifnot(max(abs(sweep(dif,1,dif[,1],'-')))<1e-7)
 cost<-acos(pmin(abs(crossprod(normcol(t(x$TR$U)%*%t(x$TR$O[[g]])),normcol(solve(w)))),1))
 map<-as.integer(solve_LSAP(cost));cost[cbind(seq_len(j$d),map)]*180/pi
 });q$angle_deg<-mean(angles)
 };q
}))
write.csv(z,file.path(base,'results.csv'),row.names=FALSE)
h<-readRDS(file.path(stage,'results/Section6_Recovery.rds'));h$angle_deg<-acos(pmin(h$cos_abs,1))*180/pi
h<-aggregate(cbind(ari,angle_deg,pi_err)~N+rep+variant,h,mean)
d<-subset(z,section=='recovery');h<-rbind(h,data.frame(N=d$N,rep=d$rep,variant='NPICA',ari=d$ARI,angle_deg=d$angle_deg,pi_err=d$pi_error));h<-aggregate(cbind(ari,angle_deg,pi_err)~N+variant,h,mean);write.csv(h,file.path(base,'recovery_comparison.csv'),row.names=FALSE);print(h)
s<-subset(z,section=='scale' & returned);h<-readRDS(file.path(stage,'results/Section6_Scale.rds'));hh<-merge(s,h,by.x=c('N','d','K','rep'),by.y=c('N','m','K','rep'));hh<-subset(hh,ok)
tt<-aggregate(cbind(wall,t_tot)~cell+method,hh,median);write.csv(tt,file.path(base,'timing_paired.csv'),row.names=FALSE);print(round(xtabs(t_tot~cell+method,tt),2));print(aggregate(wall~cell,s,median))
a<-subset(z,K==Ktrue);write.csv(do.call(rbind,lapply(split(a,paste(a$section,a$cell)),function(g)data.frame(section=g$section[1],cell=g$cell[1],n=nrow(g),returned=sum(g$returned),ARI=mean(g$ARI,na.rm=TRUE),angle=mean(g$angle_deg,na.rm=TRUE),pi_error=mean(g$pi_error,na.rm=TRUE),median_seconds=median(g$wall[g$returned]),mean_seconds=mean(g$wall[g$returned]),at_cap=sum(g$cap,na.rm=TRUE)))),file.path(base,'summary.csv'),row.names=FALSE)
