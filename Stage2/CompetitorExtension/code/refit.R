# Refit one manifest row; keep published outputs separate.
source('code/models.R')
args<-commandArgs(TRUE)
sec<-match.arg(args[1],c('section4','section5','section6'))
table<-switch(sec,section4='simulation',section5='section5',section6='section6')
jobs<-read.csv(file.path('../results/competitor-extension',paste0(table,'.csv')))
index<-as.integer(args[2]);stopifnot(length(index)==1,is.finite(index),index>=1,index<=nrow(jobs))
j<-jobs[index,];key<-switch(sec,section4=sprintf('%s_%03d',j$scenario,j$replicate),section5=sprintf('%s_%s_%03d',j$dataset,j$setting,j$replicate),section6=j$id)
d<-readRDS(file.path('replay',sec,'data',paste0(key,'.rds')))
expected<-if(sec=='section6')j$hash else j$data_hash
stopifnot(digest::digest(list(Y=d$Y,labels=d$labels),algo='sha256')==expected)
method<-if(sec=='section6')'ZH2019_ICA' else j$method
seed<-if(sec=='section5')j$seed else j$fit_seed
start<-proc.time()
a<-tryCatch(fit_new(method,d$Y,j$K,seed),error=function(e)list(error=conditionMessage(e)))
wall<-unname((proc.time()-start)[3])
ok<-is.null(a$error)&&is.finite(a$loglik)&&is.matrix(a$z)&&all(is.finite(a$z))&&all(a$z>=0)&&max(abs(rowSums(a$z)-1))<1e-6&&length(a$cluster)==nrow(d$Y)
outfile<-file.path('refits',sec,sprintf('%s_%s_K%d.rds',method,key,j$K))
atomic_rds(list(job=j,result=a,returned=ok,ARI=if(ok)mclust::adjustedRandIndex(d$labels,a$cluster) else NA_real_,wall=wall),outfile)
