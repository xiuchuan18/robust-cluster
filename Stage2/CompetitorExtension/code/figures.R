library(ggplot2)
repo<-normalizePath('../..');stage<-file.path(repo,'Stage2')
source(file.path(stage,'PlotStyle.R'))
res <- file.path(stage,'results')
out <- '../figures'

labels <- c(GMM='GMM',tMM='tMM',SNMM='MSN',STMM='MST',MSGHD='MSGHD',CGMM='CGMM',
 FW2014_MSt='MMST',ZH2019_ICA='NPICA',T2025_Manly='MMSCN',T2025_Power='PMSCN',
 MixSIT_None='MixSIT-VV',MixSIT_U='MixSIT-EV',MixSIT_Full='MixSIT-EE')
new <- names(labels)[7:10]
JOURNAL_METHOD_COLORS <- c(JOURNAL_METHOD_COLORS,'MMST'='#67539A','NPICA'='#36743E','MMSCN'='#C85B8D','PMSCN'='#836B36')
sc <- c(S1_Gaussian='S1: Gaussian',S2_Imbalanced='S2: Imbalanced',S3_GaussCopula='S3: Heavy tailed',S4_Skewed='S4: Strongly skewed',S7_CellWise10='S5: Cell-wise contamination',S6_AlignedWide='S6: Aligned, m = 10')
sets <- c(Banknote='Banknote (m = 6)',Biopsy='Biopsy (m = 9)',AIS='AIS (m = 11)',Wine='Wine (m = 13)')
src <- '../results/competitor-extension'
sim<-read.csv(file.path(src,'simulation_plot_data.csv'))
real<-read.csv(file.path(src,'section5_plot_data.csv'))
sim$Model<-factor(labels[sim$method],levels=labels);sim$Scenario<-factor(sc[sim$scenario],levels=sc)
real$Model<-factor(labels[real$method],levels=rev(labels));real$Dataset<-factor(sets[real$dataset],levels=sets)
clean<-subset(real,clean);dirty<-subset(real,!clean)
stopifnot(all(table(sim$Model,sim$Scenario)==100),all(table(clean$Dataset,clean$Model)==1),all(table(dirty$Dataset,dirty$Model)==50))
box<-function(metric,title,y,caption){
 v<-sim[[metric]];a<-sim;a$Value<-v;a<-subset(a,success & is.finite(Value))
 ggplot(a,aes(Model,Value,fill=Model))+geom_boxplot(outlier.size=.45,linewidth=.35,width=.72)+facet_wrap(~Scenario,ncol=3,scales='free_y')+scale_x_discrete(limits=if(metric=='BIC')unname(labels[labels!='NPICA']) else unname(labels))+scale_fill_methods(drop=FALSE)+labs(title=title,x=NULL,y=y,caption=caption)+theme_journal(10)+theme(axis.text.x=element_text(angle=50,hjust=1,size=8),panel.grid.major.x=element_blank(),legend.position='none',plot.caption=element_text(hjust=0))
}
style5<-function(title,note)list(facet_wrap(~Dataset,nrow=1),scale_colour_methods(drop=FALSE,guide='none'),scale_x_continuous(breaks=seq(0,1,.2)),labs(title=title,x='Adjusted Rand index',y=NULL,shape=NULL,caption=note),theme_journal(10),theme(legend.position='bottom',panel.grid.major.y=element_blank(),plot.caption=element_text(hjust=0)))
clean$status<-ifelse(clean$success & is.finite(clean$ARI),'Original data','Original fit not returned')
dirty$status<-ifelse(dirty$success & is.finite(dirty$ARI),'Contaminated fit','Fit not returned')
sentinel<-min(-.04,min(real$ARI,na.rm=TRUE)-.02)
clean$x<-ifelse(clean$status=='Original data',clean$ARI,sentinel)
dirty$x<-ifelse(dirty$status=='Contaminated fit',dirty$ARI,sentinel)
marks<-scale_shape_manual(values=c('Original data'=23,'Original fit not returned'=8,'Contaminated fit'=16,'Fit not returned'=4))
pc<-ggplot(clean,aes(x,Model,colour=Model))+geom_vline(xintercept=0,colour='grey72',linewidth=.3)+geom_point(aes(shape=status),fill='white',size=2.6,stroke=.9)+marks+expand_limits(x=c(sentinel,1))+style5('Clean real data','One fit per dataset and method. A star at the left margin indicates a missing fit.')
pd<-ggplot(dirty,aes(x,Model,colour=Model))+geom_vline(xintercept=0,colour='grey72',linewidth=.3)+geom_boxplot(data=subset(dirty,status=='Contaminated fit'),orientation='y',width=.62,outlier.shape=NA,linewidth=.45)+geom_point(aes(shape=status),alpha=.34,size=.8,position=position_jitter(height=.11,width=0,seed=2026))+geom_point(data=clean,aes(x,Model,shape=status),inherit.aes=FALSE,colour='black',fill='white',size=2,stroke=.7)+marks+style5('Cell-wise contamination','50 replicates per dataset and method. Crosses and stars at the left margin indicate missing fits.')

journal_save(file.path(out,'Figure1_ARI_Boxplots.pdf'),box('ARI',NULL,'Adjusted Rand Index',NULL),9,5.6)
journal_save(file.path(out,'Figure2_BIC_Boxplots.pdf'),box('BIC',NULL,'BIC',NULL),9,5.6)
journal_save(file.path(out,'Section5_CellWise05_ARI.pdf'),pd+labs(title=NULL,caption=NULL),9,4.5)
pl<-read.csv(file.path(res,'Section6_Scaling_Points.csv'))
raw<-subset(read.csv('../results/competitor-extension/section6.csv'),section=='scale')
stopifnot(nrow(raw)==360,sum(raw$returned)==350)
g<-subset(raw,returned)
sweeps<-list(N=subset(g,d==10 & K==3),d=subset(g,N==800 & K==3),K=subset(g,N==800 & d==10))
sn<-c(N='sample size N',d='dimension m',K='components K')
for(k in names(sweeps)){
 t<-sweeps[[k]];m<-aggregate(t$wall,list(x=t[[k]]),median)
 pl<-rbind(pl,data.frame(x=m$x,method='ZH2019_ICA',secs=m[[2]],sweep=sn[[k]],Method='NPICA'))
}
pl$Method<-factor(pl$Method,levels=c('GMM','MSGHD','NPICA','MixSIT-VV','MixSIT-EV','MixSIT-EE'))
pl$sweep<-factor(pl$sweep,levels=sn)
stopifnot(nrow(pl)==84,all(is.finite(pl$secs)),!anyNA(pl$Method))
JOURNAL_METHOD_SHAPES['NPICA']<-5
JOURNAL_METHOD_LINETYPES['NPICA']<-'longdash'
p<-ggplot(pl,aes(x,secs,colour=Method,shape=Method,linetype=Method,group=Method))+
 geom_line(linewidth=.65)+geom_point(size=1.75,stroke=.65)+facet_wrap(~sweep,scales='free_x')+
 scale_x_log10(breaks=sort(unique(pl$x)))+scale_y_log10(breaks=c(.1,.2,.5,1,2,5,10,20,50,100,200,400))+
 scale_colour_methods()+scale_shape_methods()+scale_linetype_methods()+theme_journal()+
 labs(x=NULL,y='median seconds per fit (log scale)',colour=NULL,shape=NULL,linetype=NULL)+theme(legend.position='top')
journal_save(file.path(out,'Section6_Scaling.pdf'),p,7.8,3.2)
