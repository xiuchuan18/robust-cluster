args0<-commandArgs(FALSE);this<-normalizePath(sub("^--file=","",grep("^--file=",args0,value=TRUE)[1]),winslash="/");root<-normalizePath(file.path(dirname(this),".."),winslash="/");setwd(file.path(root,"Stage2"))
source("00_CONFIG.R");source("Section4_Lib.R");source("Section5_Datasets.R");source("CalFlex_FormalLib.R")
O<-file.path(RESULTS_DIR,"calibrated-flexibility")
read_valid<-function(name){f<-file.path(O,name);if(!calflex_cache_valid(f))stop("invalid or missing: ",name);readRDS(f)}
P<-read_valid("CalFlex_PhaseMap.rds");mi<-P$method=="MSGHD"&P$ok;stopifnot(nrow(P)==8*5*50*5,length(unique(P$rep))==50,all(P$K%in%1:5),all(P$npar[mi]==mapply(calflex_msghd_npar,P$K[mi],P$d[mi],MoreArgs=list(convention="actual"))))
B<-read_valid("CalFlex_Banknote_Clean.rds");C<-read_valid("CalFlex_Banknote_CellWise05.rds");stopifnot(nrow(B)==9,nrow(C)==50*9,length(unique(C$rep))==50,all(C$alpha==.05),all(C$eta==5))
Z<-read_valid("CalFlex_BayesCeiling.rds");stopifnot(nrow(Z)==6,all(Z$n_rep==100),all(is.finite(Z$ARI_Bayes)))
for(f in c("CalFlex_PhaseMap_Summary.csv","CalFlex_PhaseMap_Paths.csv","CalFlex_Banknote_Summary.csv","CalFlex_CeilingAttainment_Profile.csv"))if(!file.exists(file.path(O,f)))stop("missing summary: ",f)
cat("calibrated-flexibility integrity audit passed\n")
