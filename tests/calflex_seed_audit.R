args0<-commandArgs(FALSE);this<-normalizePath(sub("^--file=","",grep("^--file=",args0,value=TRUE)[1]),winslash="/");root<-normalizePath(file.path(dirname(this),".."),winslash="/");setwd(file.path(root,"Stage2"))
source("00_CONFIG.R");source("Section4_Lib.R");source("Section5_Datasets.R");source("CalFlex_FormalLib.R")
cells<-data.frame(mechanism=c("cellSD","cellSD","cellSD","cellSD","cellSD","cellSD","tailMM","skewMM"),severity=c(0,2/3,1,0,2/3,1,1,1),separation=c(4,4,4,2.5,2.5,2.5,2.5,2.5))
keys<-character();seeds<-integer();add<-function(key,seed){keys<<-c(keys,key);seeds<<-c(seeds,seed)}
for(i in seq_len(nrow(cells)))for(r in 1:50){j<-cells[i,];add(paste("phase-data",j$mechanism,j$severity,j$separation,r),calflex_seed("PhaseMap",j$mechanism,j$severity,j$separation,"data",r));for(m in calflex_phase_methods)for(K in 1:5){tag<-paste(CALFLEX_SEED_PREFIX,"PhaseMap",j$mechanism,j$severity,j$separation,sep="/");rp<-paste(m,r,"K",K,sep="_");add(paste("phase-fit",i,m,r,K),calflex_fit_seed(tag,m,rp))}}
for(m in calflex_methods){tag<-paste(CALFLEX_SEED_PREFIX,"Banknote","clean",sep="/");add(paste("bank-clean",m),calflex_fit_seed(tag,m,1L))}
for(r in 1:50){add(paste("bank-contam",r),calflex_seed("Banknote","cellwise",.05,5,"contam",r));for(m in calflex_methods){tag<-paste(CALFLEX_SEED_PREFIX,"Banknote","cellwise05_eta5",sep="/");add(paste("bank-fit",m,r),calflex_fit_seed(tag,m,r))}}
if(anyNA(seeds)||any(seeds<=0))stop("invalid formal seed")
dup<-duplicated(seeds)|duplicated(seeds,fromLast=TRUE)
if(any(dup)){print(data.frame(key=keys[dup],seed=seeds[dup]));stop("formal seed collision")}
cat(sprintf("calibrated-flexibility seed audit passed: %d unique semantic streams\n",length(seeds)))
