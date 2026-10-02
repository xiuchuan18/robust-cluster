# Recompute the exact six-cell Bayes ceilings using the production seeds, then
# combine them with the current formal 100-replicate Section 4 fits.
setwd(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/")))
source("00_CONFIG.R");source("Section4_Lib.R");source("Section5_Datasets.R");source("CalFlex_FormalLib.R")
OUT<-file.path(RESULTS_DIR,"calibrated-flexibility");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
formal_file<-file.path(OUT,"CalFlex_BayesCeiling.rds")
raw_file<-file.path(RESULTS_DIR,"Section4_BayesCeiling.rds")
raw<-if(file.exists(raw_file))tryCatch(readRDS(raw_file),error=function(e)NULL)else NULL
raw_p<-if(!is.null(raw))attr(raw,"provenance",exact=TRUE)else NULL
fresh_raw<-is.data.frame(raw)&&nrow(raw)==6L&&all(raw$n_rep==100L)&&
  is.list(raw_p)&&identical(raw_p$run_id,EXPERIMENT_RUN_ID)&&
  identical(raw_p$scheme,RNG_SCHEME_VERSION)&&
  identical(raw_p$r_version,"R version 4.4.1 (2024-06-14 ucrt)")
if(!calflex_cache_valid(formal_file)&&!fresh_raw){BAYES_N_REP<-100L;source("Section4_BayesCeiling.R")}else if(fresh_raw)cat("reusing current completed Bayes recomputation\n")
suppressPackageStartupMessages(library(dplyr))
ceil<-if(calflex_cache_valid(formal_file))readRDS(formal_file)else readRDS(raw_file)
sim_file<-file.path(RESULTS_DIR,"Section4_Simulation_Results.rds")
if(!cache_valid(sim_file))stop("current formal Section 4 result is absent or stale: ",sim_file)
sim<-readRDS(sim_file)
long<-dplyr::bind_rows(lapply(names(sim),function(cell){z<-sim[[cell]];ac<-grep("_ARI$",names(z),value=TRUE);dplyr::bind_rows(lapply(ac,function(v)data.frame(cell=cell,rep=z$sim_id,method=sub("_ARI$","",v),ARI=z[[v]])))}))%>%dplyr::left_join(ceil%>%dplyr::select(cell,ARI_Bayes),by="cell")%>%dplyr::mutate(Method=unname(calflex_method_label[method]),attainment=ARI/ARI_Bayes)
cell<-long%>%group_by(Method,cell)%>%summarise(n_returned=sum(is.finite(ARI)),mean_ARI=mean(ARI,na.rm=TRUE),mean_attainment=mean(attainment,na.rm=TRUE),sd_attainment=sd(attainment,na.rm=TRUE),.groups="drop")
profile<-cell%>%group_by(Method)%>%summarise(cells=n(),average_attainment=mean(mean_attainment),worst_cell_attainment=min(mean_attainment),attainment_sd_across_cells=sd(mean_attainment),total_returned=sum(n_returned),.groups="drop")%>%arrange(desc(worst_cell_attainment))
calflex_save(ceil,file.path(OUT,"CalFlex_BayesCeiling.rds"));calflex_write_csv(ceil,file.path(OUT,"CalFlex_BayesCeiling.csv"));calflex_write_csv(cell,file.path(OUT,"CalFlex_CeilingAttainment_ByCell.csv"));calflex_write_csv(profile,file.path(OUT,"CalFlex_CeilingAttainment_Profile.csv"));cat("Bayes ceiling/profile complete\n")
