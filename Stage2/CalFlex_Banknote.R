# Banknote negative-control dataset under the production clean and semi-real
# protocols.  The semi-real arm uses alpha=.05, eta=5, and 50 paired copies.
args0<-commandArgs(FALSE);SCRIPT<-normalizePath(sub("^--file=","",grep("^--file=",args0,value=TRUE)[1]),winslash="/");setwd(dirname(SCRIPT))
suppressWarnings(suppressMessages({source("00_CONFIG.R");source_stage1();source("Section4_Lib.R");source("Section5_Datasets.R");source("CalFlex_FormalLib.R")}))
OUT<-file.path(RESULTS_DIR,"calibrated-flexibility");CACHE<-file.path(OUT,"cache","banknote");LOG<-file.path(OUT,"logs","banknote")
dir.create(CACHE,recursive=TRUE,showWarnings=FALSE);dir.create(LOG,recursive=TRUE,showWarnings=FALSE)
ALPHA<-.05;ETA<-5;NREP<-50L
cache_file<-function(rep)file.path(CACHE,sprintf("banknote_contam_r%02d.rds",rep))

fit_one<-function(method,X,d,tag,rep){
 fit_tag<-paste(CALFLEX_SEED_PREFIX,"Banknote",tag,sep="/");seed<-calflex_fit_seed(fit_tag,method,rep);set.seed(seed);t0<-proc.time()[3]
 rr<-tryCatch(fit_multistart(method,X,d$K,R_RESTART,COMPETITOR_MAX_ITER,
   tag=fit_tag,rep=rep),error=function(e)NULL)
 sec<-proc.time()[3]-t0;ok<-!is.null(rr)&&!is.null(rr$fit)&&!is.null(rr$cluster)
 ic<-if(ok)calc_information_criteria(X,rr$fit,method,d$K,d$m)else list(AIC=NA_real_,BIC=NA_real_,ICL=NA_real_)
 data.frame(method,success=ok,converged=calflex_converged(rr,method),
  ARI=if(ok)calc_ari(d$y,rr$cluster)else NA_real_,NMI=if(ok)calc_nmi(d$y,rr$cluster)else NA_real_,
  AIC=ic$AIC,BIC=ic$BIC,ICL=ic$ICL,seconds=sec,fit_seed=seed)
}

run_clean<-function(){d<-calflex_load_data("Banknote");z<-do.call(rbind,lapply(calflex_methods,function(m)fit_one(m,d$X,d,"clean",1L)));cbind(dataset="Banknote",setting="clean",rep=1L,n=d$n,d=d$m,K=d$K,z)}
run_contam<-function(rep){
 d<-calflex_load_data("Banknote");seed<-calflex_seed("Banknote","cellwise",ALPHA,ETA,"contam",rep);set.seed(seed)
 mask<-matrix(runif(length(d$X))<ALPHA,nrow(d$X),ncol(d$X));shock<-matrix(rnorm(length(d$X),sd=ETA),nrow(d$X),ncol(d$X));Xc<-d$X+mask*shock
 z<-do.call(rbind,lapply(calflex_methods,function(m)fit_one(m,Xc,d,"cellwise05_eta5",rep)))
 cbind(dataset="Banknote",setting="cellwise05_eta5",rep=rep,alpha=ALPHA,eta=ETA,
  contamination_seed=seed,n=d$n,d=d$m,K=d$K,n_contaminated=sum(mask),
  fraction_contaminated=mean(mask),n_rows_affected=sum(rowSums(mask)>0),
  fraction_rows_affected=mean(rowSums(mask)>0),z)
}

args<-commandArgs(TRUE)
if(length(args)==2L&&args[1]=="rep"){r<-as.integer(args[2]);f<-cache_file(r);if(!calflex_cache_valid(f))calflex_save(run_contam(r),f);quit(save="no")}
if(!length(args)||args[1]!="driver")stop("use rep <r> or driver [parallel]")
np<-if(length(args)>=2)as.integer(args[2])else 8L
clean_file<-file.path(OUT,"CalFlex_Banknote_Clean.rds");if(!calflex_cache_valid(clean_file))calflex_save(run_clean(),clean_file)
reps<-seq_len(NREP);todo<-reps[!vapply(lapply(reps,cache_file),calflex_cache_valid,logical(1))]
cat(sprintf("formal Banknote contamination: %d/%d reps outstanding; %d concurrent\n",length(todo),NREP,np))
if(length(todo)){active<-list();nx<-1L;launch<-function(i){r<-todo[i];processx::process$new(file.path(R.home("bin"),"Rscript.exe"),c(SCRIPT,"rep",r),stdout=file.path(LOG,sprintf("contam_r%02d.log",r)),stderr="2>&1",cleanup=TRUE,env=c(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))}
 while(nx<=length(todo)||length(active)){while(nx<=length(todo)&&length(active)<np){active[[as.character(nx)]]<-launch(nx);nx<-nx+1L};Sys.sleep(.5);fin<-names(active)[!vapply(active,function(p)p$is_alive(),logical(1))];for(id in fin){p<-active[[id]];if(p$get_exit_status()!=0L)warning("Banknote job failed: ",id);active[[id]]<-NULL}}}
files<-lapply(reps,cache_file);if(!all(vapply(files,calflex_cache_valid,logical(1))))stop("Banknote contamination incomplete")
D<-do.call(rbind,lapply(files,readRDS));calflex_save(D,file.path(OUT,"CalFlex_Banknote_CellWise05.rds"));calflex_write_csv(readRDS(clean_file),file.path(OUT,"CalFlex_Banknote_Clean.csv"));calflex_write_csv(D,file.path(OUT,"CalFlex_Banknote_CellWise05.csv"));cat("formal Banknote experiment complete\n")
