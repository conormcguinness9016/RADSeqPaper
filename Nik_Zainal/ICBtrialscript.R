library(XML)
library(RCurl)
library(stringr)
library(devtools)
library(stringr)
library(tidyverse)
library(SomaticSignatures)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(reshape2)
library(bedr)
library(GenomicRanges)
library(data.table)
library(cowplot)
#library(emdbook)
library(data.table)
library(parallel)
library(ggsignif)
library(gridExtra)
library(gtable)
library(grid)
library(scales)
library(pROC)
library(ggpubr)
AllEnzInfo<-readRDS("AllEnzInfo.RDS")
libs<-c("F1","SpeI", "KpnI")
AllEnzInfofilt<-AllEnzInfo[names(AllEnzInfo) %in% libs]
Ahn<-length(which(AllEnzInfofilt[[1]]$APOBECHyper=="APOBEC_hyper"))
Ahn
set.seed(10)
random_ss<-round(runif(Ahn,min=1,max=nrow(AllEnzInfofilt[[1]])))
random_ss[duplicated(random_ss)]<-random_ss[duplicated(random_ss)]+10

midpoint<-c(12,8,10)
minresp<-c(0.04,0.04,0.04)
maxresp<-c(0.5,0.5,0.5)
sd=c(8,2,0)
model=c("Gradual","Moderate","Extreme")
dir.create("ICBtrials")
model_df<-data.frame(midpoint,minresp,maxresp,sd=sd,model=model)

ResponsePlot<-list()
ORplot<-list()
prop_plot<-list()
proptrialLibs<-list()
ORplotTrialLibs<-list()
m<-1
while(m<nrow(model_df)+1){
  mod<-model_df[m,]
  midpoint<-mod$midpoint
  minresp<-mod$minresp
  maxresp<-mod$maxresp
  sd=mod$sd
  model=mod$model
  folder<-paste0("ICBtrials/",model,"/")
  print(folder)
  dir.create(folder)
  dist<-seq(min(AllEnzInfofilt[[1]]$mutratewgs), max(AllEnzInfofilt[[1]]$mutratewgs)+1, by =0.01)
  dfdist<-data.frame(mutratewgs=dist,prob=rescale(pnorm(dist,midpoint,sd),c(minresp,maxresp)))
  ggplot(dfdist,aes(mutratewgs, prob))+geom_point()+ylim(c(0,1))+xlim(c(0,30))+
  xlab("Mutation rate (mutations/Mb)")+ylab("Probability of \nresponse to treatment")+
  theme(axis.text=element_text( size =20), axis.text.x = element_text(hjust = 1,vjust=0.5),axis.title = element_text(size=25))
  ggsave(paste0(folder,"modelplot.png"))

cutofftests<-function(x,df){
  cutoff<-x
  df$WGSpred[df$mutratewgs>cutoff]<-1
  df$WGSpred[df$mutratewgs<cutoff]<-0
  df$GBSpred[df$mutrategbs>cutoff]<-1
  df$GBSpred[df$mutrategbs<cutoff]<-0
  df$WGSpredchar[df$mutratewgs>cutoff]<-"Response"
  df$WGSpredchar[df$mutratewgs<cutoff]<-"Non-response"
  df$GBSpredchar[df$mutrategbs>cutoff]<-"Response"
  df$GBSpredchar[df$mutrategbs<cutoff]<-"Non-response"
  df$actual<-factor(df$actual, levels=c("Non-response","Response"))
  df$WGSpred<-factor(df$WGSpred, levels=0:1)
  df$GBSpred<-factor(df$GBSpred,levels=0:1)
  gbstab<-table(actual=df$actual,df$GBSpred)
  gbstab<-table(actual=df$actual,df$GBSpred)
  gbsproptab<-prop.table(gbstab)
  fishgbs<-fisher.test(gbstab)
  gbs_pval<-fishgbs$p.value
  gbs_OR<-fishgbs$estimate
  gbs_OR_lwr<-fishgbs$conf.int[1]
  gbs_OR_upr<-fishgbs$conf.int[2]
  PPVgbs<-gbstab[2,2]/sum(gbstab[,2])
  NPVgbs<-gbstab[1,1]/sum(gbstab[,1])
  specificitygbs<-(gbstab[1,1]/sum(gbstab[1,]))
  sensitivitygbs<-gbstab[2,2]/sum(gbstab[2,])
  df$GBSpred<-as.numeric(df$GBSpred)
  if(length(unique(df$actual))>1){
    rocObj_gbs=roc(df$actual,df$GBSpred,ci=T)
    AUC=as.numeric(rocObj_gbs$auc)}else{
      AUC=0.5
    }
  gbsdf<-data.frame(TN=gbsproptab[1,1],FN=gbsproptab[2,1],
                    TP=gbsproptab[2,2],FP=gbsproptab[1,2],
                    PPV=PPVgbs, NPV=NPVgbs,Spec=specificitygbs,Sens=sensitivitygbs,
                    AUC=AUC,
                    OR=gbs_OR,pval=gbs_pval,OR_lwr=gbs_OR_lwr, OR_upper=gbs_OR_upr,
                    cutoff = x, Library=df$Library[1])
  gbsdf
}


LibraryStatsICBtrials<-function(df,midpoint,sd,minresp,maxresp, lib,number_reps, subset=NULL){
  df<-df[names(df) %in% lib]
  dfWGS<-df[[1]]
  dfWGS$mutrategbs<-dfWGS$mutratewgs
  dfWGS$Library<-"WGS"
  prob_wgs<-rescale(pnorm(dfWGS$mutratewgs,midpoint,sd),c(minresp,maxresp))
  simresp<-rbinom(prob_wgs,1,prob_wgs)
  df<-c(list(dfWGS),df)
  dfAh<-list()
  dfAll<-list()
  valsAll<-list()
  valsAh<-list()
  vals<-list()
  sumdf<-list()
  for(i in 1:length(df)){
    df[[i]]$prob_wgs<-prob_wgs
    df[[i]]$simresp<-simresp
    df[[i]]$prob_gbs<-rescale(pnorm(df[[i]]$mutrategbs,midpoint,sd),c(minresp,maxresp))
    df[[i]]$actual[df[[i]]$simresp ==1]<-"Response"
    df[[i]]$actual[df[[i]]$simresp ==0]<-"Non-response"
    dfAh[[i]]<-subset(df[[i]], APOBECHyper=="APOBEC_hyper")
    dfAh[[i]]$Cohort<-"APOBEC_hyper"
    dfAll[[i]]<-df[[i]][random_ss,]
    dfAll[[i]]$Cohort<-"Unselected"
    valsAll[[i]]<-bind_rows(lapply(seq(from=1,to=15),function(x){
      cutofftests(x,dfAll[[i]])}))
    valsAh[[i]]<-bind_rows(lapply(seq(from=1,to=15),function(x){
      cutofftests(x,dfAh[[i]])}))
    valsAll[[i]]$pval.adj<-p.adjust(valsAll[[i]]$pval)
    valsAh[[i]]$pval.adj<-p.adjust(valsAh[[i]]$pval)
    valsAh[[i]]$Cohort<-"APOBEC_hyper"
    valsAll[[i]]$Cohort<-"Unselected"
    vals[[i]]<-bind_rows(valsAh[[i]],valsAll[[i]])
    sumdf[[i]]<-bind_rows(dfAh[[i]],dfAll[[i]])
  }
  sumdf<-bind_rows(sumdf)
  vals<-bind_rows(vals)
  return(list(sumdf,vals))
  }


r=1
all_stats<-data.frame()
all_profiles<-data.frame()
while(r<100){
  print(paste("r=",r))
  set.seed(r*1000)
  lsICB_extreme<-LibraryStatsICBtrials(df=AllEnzInfofilt,
                                       midpoint=midpoint, 
                                       sd=sd, 
                                       minresp=minresp,
                                       maxresp=maxresp,
                                       lib=libs,
                                       number_reps=1)
  profiles_extreme<-lsICB_extreme[[1]]
  profiles_extreme$trial=r
  stats_extreme<-lsICB_extreme[[2]]
  stats_extreme$trial=r
  all_profiles<-bind_rows(all_profiles,profiles_extreme)
  all_stats<-bind_rows(all_stats,stats_extreme)
  r<-r+1
}
profiles_extreme<-all_profiles
prop.table(table(Response=profiles_extreme$simresp,Cohort=profiles_extreme$Cohort),2)
stats_extreme<-all_stats
stats_extreme$sig[stats_extreme$pval.adj<0.05]<-"Significant"
stats_extreme$sig[stats_extreme$pval.adj>0.05]<-"Not Significant"
stats_extreme
ResponsePlot[[m]]<-ggplot(profiles_extreme,aes(mutratewgs, prob_wgs))+
  geom_point(data=dfdist,aes(mutratewgs, prob),size=0.2,col="black")+ylim(c(0,1))+
  geom_jitter(aes(col=actual),alpha=0.25,width=0.1,size=3)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  guides(col=guide_legend(title="Response"))+
  xlab("Mutation rate (mutations/Mb)")+ylab("Probability of \nresponse to treatment")
ResponsePlot[[m]]
ggsave(paste0(folder,"ResponsePlot.png"),width=15,height=10)
ggplot(profiles_extreme,aes(prob_gbs, prob_wgs))+geom_point(aes(col=actual),alpha=0.5)+ylim(c(0,1))+xlim(c(0,1))+
  geom_abline(intercept=0,slope=1)
saveRDS(stats_extreme, paste0(folder,"stats_extreme.RDS"))
saveRDS(profiles_extreme, paste0(folder,"profiles_extreme.RDS"))
stats_extreme<-readRDS(paste0(folder,"stats_extreme.RDS"))
profiles_extreme<-readRDS(paste0(folder,"profiles_extreme.RDS"))
stats_extremeWGS<-subset(stats_extreme, Library=="WGS")
sigdf<-NULL
for(i in 1:15){
  tdf<-subset(stats_extremeWGS, cutoff ==i)
  sum<-summarise(group_by(tdf, Cohort), sd(AUC),sd(OR),)
  sdAh=as.numeric(sum[1,2])
  sdUns=as.numeric(sum[2,2])
  sigtab<-table(factor(tdf$sig,levels=c("Not Significant", "Significant")),tdf$Cohort)
  ch<-chisq.test(sigtab)
  propAh_sig<-prop.table(sigtab,2)[2,1]
  propUns_sig<-prop.table(sigtab,2)[2,2]
  t<-t.test(AUC ~ Cohort,tdf, paired=TRUE)
  sigdf[[i]]<-data.frame(t=t$parameter,
                         pval=t$p.value, 
                         meanAh=t$estimate[1],
                         meanUns=t$estimate[2],
                         sdAh=sdAh,
                         sdUns=sdUns,
                         cutoff=i,chisquare=ch$statistic,chisquarep=ch$p.value,
                         propAh_sig=propAh_sig,propUns_sig=propUns_sig)
}
sigdf<-bind_rows(sigdf)
sigdf$pval.adj<-p.adjust(sigdf$pval)
sigdf$chisquarepadj<-p.adjust(sigdf$chisquarep)
saveRDS(sigdf, paste0(folder,"sigdf.RDS"))
# ggplot(stats_extremeWGS,aes(as.factor(cutoff),PPV,col=Cohort))+ylim(c(0,1))+geom_boxplot(outlier.shape=NA)+
#   geom_point(position = position_dodge(width=0.75),alpha=0.1)+ylim(c(0,1))+
#   scale_colour_manual(values=c("red", "black"))+
#   #ylim(-5,5)+
#   theme(axis.text=element_text(angle=90, size =10), axis.text.x = element_text(hjust = 1,vjust=0.5),axis.title = element_text(size=15))+
#   ylab("PPV for cutoff")+xlab("TMB Cutoff")
ggplot(stats_extremeWGS, aes(as.factor(cutoff), AUC,col=Cohort))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width=0.75),alpha=0.1)+ylim(c(0,1))+
  scale_colour_manual(values=c("red", "black"))+
  #ylim(-5,5)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("AUC for cutoff")+xlab("TMB Cutoff")
ggsave(paste0(folder,"AUCplotCohorts.png"))


# ggplot(stats_extremeWGS, aes(Spec, (1-Sens),col=Cohort))+geom_point()+ylim(c(0,1))
stats_extremeWGS$OR[!is.finite(stats_extremeWGS$OR)]<-max(stats_extremeWGS$OR[is.finite(stats_extremeWGS$OR)])
ORplot[[m]]<-ggplot(stats_extremeWGS, aes(as.factor(cutoff),OR,col=Cohort))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width=0.75),alpha=0.1)+
  scale_colour_manual(values=c("red", "black"))+
  #ylim(-5,5)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("OR for response vs \nnon-response at cutoff")+xlab("TMB Cutoff")
ORplot[[m]]
ggsave(paste0(folder,"ORpvalplotCohorts.png"))
ggplot(stats_extremeWGS, aes(as.factor(cutoff), -log10(pval.adj), col =Cohort))+geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept=-log10(0.05)))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width=0.75),alpha=0.1)+
  scale_colour_manual(values=c("red", "black"))+
  #ylim(-5,5)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("Adjusted -log10p for OR at cutoff")+xlab("TMB Cutoff")
ggsave(paste0(folder,"ORsignificanceplotCohorts.png"))
propdf<-data.frame(Cohort=rep(c("APOBEC_hyper","Unselected"), each=nrow(sigdf)), proportionsig=c(sigdf$propAh_sig,sigdf$propUns_sig), cutoff=rep(sigdf$cutoff,2))
prop_plot[[m]]<-ggplot(propdf, aes(factor(cutoff),proportionsig, fill=Cohort))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=c("red", "black"))+
  geom_text(data=sigdf,aes(factor(cutoff), label=signif(chisquarepadj,3), y=0.75,angle=90, size =5), inherit.aes = FALSE, size=5)+
  ylim(0,1)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("Proportion of \ntrials reporting \nsignificant association")+xlab("TMB Cutoff")
prop_plot[[m]]
ggsave(paste0(folder,"ORsignificanceplotCohorts.png"))


cohorts<-unique(stats_extreme$Cohort)
dir.create(paste0(folder,"cohorts"))
c<-1
while(c<(length(cohorts)+1)){
  dir.create(paste0(folder,"cohorts/",cohorts[c]))
  stats_extremeCohort<-subset(stats_extreme, Cohort==cohorts[c])
  stats_extremeCohort$misclass_rate<-stats_extremeCohort$FN+stats_extremeCohort$FP
ggplot(stats_extremeCohort, aes(as.factor(cutoff), AUC,col=Library))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width=0.75),alpha=0.1)+ylim(c(0,1))+
  #scale_colour_manual(values=c("red", "black"))+
  #ylim(-5,5)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("AUC for cutoff")+xlab("TMB Cutoff")
ggsave(paste0(folder,"cohorts/",cohorts[c],"/AUCplotLibraries.png"))

stats_extremeCohort$OR[!is.finite(stats_extremeCohort$OR)]<-max(stats_extremeCohort$OR[is.finite(stats_extremeCohort$OR)])
ORplotTrialLibs[[c]]<-ggplot(stats_extremeCohort, aes(as.factor(cutoff),OR,col=Library))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width=0.75),alpha=0.1)+
  #scale_colour_manual(values=c("red", "black"))+
  #ylim(-5,5)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("OR for response vs \nnon-response at cutoff")+xlab("TMB Cutoff")
ORplotTrialLibs[[c]]
ggsave(paste0(folder,"cohorts/",cohorts[c],"/ORplotLibraries.png"))
ggplot(stats_extremeCohort, aes(as.factor(cutoff), -log10(pval.adj), col =Library))+geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept=-log10(0.05)))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width=0.75),alpha=0.1)+
  #scale_colour_manual(values=c("red", "black"))+
  #ylim(-5,5)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("Adjusted -log10p for OR at cutoff")+xlab("TMB Cutoff")
ggsave(paste0(folder,"cohorts/",cohorts[c],"/pvalplotLibraries.png"))
sig_df<-NULL
fish<-NULL
for(i in 1:max(stats_extremeCohort$cutoff)){
  print(i)
  ssA<-subset(stats_extremeCohort, cutoff ==i)
  sigtab<-table(factor(ssA$sig, levels=c("Significant", "Not Significant")), factor(ssA$Library))
  fish[[i]]<-fisher.test(sigtab)
  sigdf<-data.frame(prop.table(sigtab,2), cutoff =i,fisherpval=fish[[i]]$p.value)
  sig_df[[i]]<-sigdf[sigdf$Var1=="Significant",]}
sig_df<-bind_rows(sig_df)
pval.adjust<-p.adjust(sig_df$fisherpval[!duplicated(sig_df$cutoff)])
chsqdf<-data.frame(cutoff=1:15,pval.adjust)
proptrialLibs[[c]]<-ggplot(sig_df, aes(factor(cutoff),Freq, fill=Var2))+geom_bar(stat="identity",position="dodge")+#scale_colour_manual(values=c("red", "black"))+
  ylim(0,1)+
  theme(axis.text=element_text(size =15), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_text(size=18))+
  ylab("Proportion of \ntrials reporting \nsignificant association")+xlab("TMB Cutoff")+
  geom_text(data=chsqdf,aes(factor(cutoff),label=round(pval.adjust,3),y=0.75),angle=90,inherit.aes = FALSE,size=5)+
  guides(fill=guide_legend(title="Library"))
proptrialLibs[[c]]
ggsave(paste0(folder,"cohorts/",cohorts[c],"/proptrialsLibraries.png"))
c<-c+1
}

p<-plot_grid(ORplot[[m]]+theme(legend.position = "none"),prop_plot[[m]])
p
ggsave(plot=p,paste0(folder,"SummaryPlot.png"))
g<-plot_grid(ORplotTrialLibs[[1]]+theme(legend.position = "none"),proptrialLibs[[1]])
TrialPlot<-ggarrange(ORplot[[m]]+theme(legend.position = "none"),prop_plot[[m]]+theme(legend.position = "none"), 
          ORplotTrialLibs[[1]]+theme(legend.position = "none"),proptrialLibs[[1]]+theme(legend.position = "none"),ncol=2,nrow=2)
TrialPlotResponse<-ggarrange(ResponsePlot[[m]]+theme(legend.position="none"),TrialPlot,heights=c(1,2),ncol=1)
TrialPlotResponse
ggsave(paste0(folder,"TrialSummaryPlot.png"),width=15,height=12)
ggsave(plot=p,paste0(folder,"CohortSummaryPlot.png"),width=15,height=10)
ggsave(plot=g,paste0(folder,"LibrarySummaryPlot.png"),width=15,height=10)
m<-m+1
}


# df<-subset(AllEnzInfobound, Library=="F1")
# dfAh<-subset(AllEnzInfobound, Library=="F1" & APOBECHyper=="APOBEC_hyper")
# cutofftests_misclass<-function(x,df){
#   cutoff<-x
#   df$WGSpred[df$mutratewgs>cutoff]<-1
#   df$WGSpred[df$mutratewgs<cutoff]<-0
#   df$GBSpred[df$mutrategbs>cutoff]<-1
#   df$GBSpred[df$mutrategbs<cutoff]<-0
#   df$WGSpredchar[df$mutratewgs>cutoff]<-"Treat"
#   df$WGSpredchar[df$mutratewgs<cutoff]<-"Not treat"
#   df$GBSpredchar[df$mutrategbs>cutoff]<-"Treat"
#   df$GBSpredchar[df$mutrategbs<cutoff]<-"Not treat"
#   misclasstab<-table(WGS_predict=df$WGSpredchar,GBS_predict=df$GBSpredchar)
#   PPV<-misclasstab[2,2]/sum(misclasstab[,2])
#   NPV<-misclasstab[1,1]/sum(misclasstab[,1])
#   Acc<-(misclasstab[1,1]+misclasstab[2,2])/sum(misclasstab)
#   Sens<-misclasstab[2,2]/sum(misclasstab[2,])
#   Spec<-misclasstab[1,1]/sum(misclasstab[1,])
#   misclassrate<-1-Acc
#   roc<-roc(df$WGSpred, df$GBSpred)
#   AUC<-as.numeric(roc$auc)
#   fishgbs<-fisher.test(misclasstab)
#   data.frame(cutoff=x, Acc,PPV,
#              NPV,Sens,Spec,AUC,pfish=fishgbs$p.value,misclassrate)
# }
# misclassdf<-bind_rows(lapply(seq(15), function(x){
#   dfUns<-cutofftests_misclass(x,df=df)
#   dfUns$Cohort<-"Unselected"
#   dfUns
#   dfAPO<-cutofftests_misclass(x,df=dfAh)
#   dfAPO$Cohort<-"APOBEC"
#   bind_rows(dfAPO,dfUns)
# }))
# ggplot(misclassdf, aes(cutoff,Acc,fill=Cohort))+geom_bar(stat="identity",position="dodge")+ylim(0,1)
# ggplot(misclassdf, aes(cutoff,PPV,fill=Cohort))+geom_bar(stat="identity",position="dodge")+ylim(0,1)
# ggplot(misclassdf, aes(cutoff,Sens,fill=Cohort))+geom_bar(stat="identity",position="dodge")+ylim(0,1)
# ggplot(misclassdf, aes(cutoff,Spec,fill=Cohort))+geom_bar(stat="identity",position="dodge")+ylim(0,1)
# ggplot(misclassdf, aes(cutoff,NPV,fill=Cohort))+geom_bar(stat="identity",position="dodge")+ylim(0,1)
# ggplot(misclassdf, aes(cutoff, AUC, fill = Cohort))+geom_bar(stat="identity",position = "dodge")+ylim(0,1)
# ggplot(misclassdf, aes(cutoff, misclassrate,fill=Cohort))+geom_bar(stat="identity",position = "dodge")+ylim(0,1)
