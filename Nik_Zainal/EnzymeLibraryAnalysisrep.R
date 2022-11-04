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
# REfraglibs<-readRDS("REfraglibs.RDS")
AllEnzInfo<-readRDS("AllEnzInfofilt.RDS")
dfprops<-readRDS("dfprops.RDS")
covs<-readRDS("covs.RDS")
takeoff<-readRDS("takeoff.RDS")
MLcors<-readRDS("MLcors.RDS")
wgscov<-sum(seqlengths(Hsapiens)[1:24])

enz<-names(AllEnzInfo)
covs<-covs[covs$Library %in% enz,]
enz
reps<-covs$Library[round(seq(1,nrow(covs),length.out = 9))]
reps[9:12]<-c("MspI","F1","KpnI","SpeI")
reps
takeoff<-takeoff[takeoff$Library %in% reps,]
AllEnzInfobound<-bind_rows(AllEnzInfo)
AllEnzInfobound<-AllEnzInfobound[AllEnzInfobound$Library %in% reps,]
MLcors<-MLcors[MLcors$Library %in% reps,]
takeoff_All<-subset(takeoff, APOBEChyper == "All_samples")
prop_plot_enzymes_all<-ggplot(takeoff_All, aes(reorder(Library, mean_nmutsgbs), prop)) + geom_bar(stat="identity", position="dodge")+#geom_point()+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  #geom_bar(stat="identity", position ="dodge")+
  # geom_text(aes(x=first_sig_Ah,y=1.02, colour =APOBEChyper[1], label="*"), size =10, show.legend = FALSE)+
  geom_text(aes(x=reorder(Library, mean_nmutsgbs),y=1.02, label=sig), size =25,show.legend = FALSE)+
  xlab("Library")+
  ylab("Proportion of samples \nresulting in >90% \ncosine similarity") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 25), axis.text =  element_text(size = 22),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 22), legend.title = element_text(size = 25))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"))
prop_plot_enzymes_all
dir.create("figs")
dir.create("figs/Enzymesrep/")
#ggsave(plot=prop_plot_enzymes_all,"figs/Enzymesrep/Proportion_plot_All_Enzymes.png",width=20,height=8)

takeoff_APO<-subset(takeoff, APOBEChyper != "All_samples")
prop_plot_enzymes_APO<-ggplot(takeoff_APO, aes(reorder(Library, mean_nmutsgbs), prop,fill=APOBEChyper)) + geom_bar(stat="identity", position="dodge")+#geom_point()+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  #geom_bar(stat="identity", position ="dodge")+
  # geom_text(aes(x=first_sig_Ah,y=1.02, colour =APOBEChyper[1], label="*"), size =10, show.legend = FALSE)+
  geom_text(aes(x=reorder(Library, mean_nmutsgbs),y=1.02, label=sig, col =APOBEChyper), size =25,show.legend = FALSE)+
  xlab("Library")+
  ylab("Proportion of samples \nresulting in >90% \ncosine similarity") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 25), axis.text =  element_text(size = 22),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 22), legend.title = element_text(size = 25))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"))
prop_plot_enzymes_APO
#ggsave(plot=prop_plot_enzymes_APO,"figs/Enzymesrep/Proportion_plot_Apo_Enzymes.png",width=22,height=8)


AllEnzCosineplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), cosimdiag)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = 0.2,width=0.08,height=0,size=5)+
  scale_colour_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 25), axis.text =  element_text(size = 22),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 22), legend.title = element_text(size = 25))+
  ylab("Cosine similarity \n to WGS profile")+xlab("Library")
#ggsave(plot=AllEnzCosineplot, "figs/Enzymesrep/AllEnzCosineplot.png",width=22,height=8)
AllEnzCosineplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), cosimdiag,col=APOBECHyper)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.08),alpha = 0.2)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 25), axis.text =  element_text(size = 22),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 22), legend.title = element_text(size = 25))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"),
         override.aes=list(alpha=1))+ylab("Cosine similarity \n to WGS profile")+xlab("Library")
AllEnzMutsplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), nmutsgbs)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = 0.25,width=0.1)+
  scale_colour_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 25), axis.text =  element_text(size = 22),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 22), legend.title = element_text(size = 25))+
  ylab("Number of mutations \ndetected (log10 scale)")+xlab("Library")+
  scale_y_continuous(trans='log10')
#ggsave(plot=AllEnzMutsplot, "figs/Enzymesrep/AllEnzMutsplot.png",width=22, height =8)
#ggsave(plot=AllEnzCosineplotAPO, "figs/Enzymesrep/AllEnzCosineplotAPO.png",width=22,height=8)
Fig3<-plot_grid(AllEnzMutsplot,AllEnzCosineplot,prop_plot_enzymes_all,align="v", axis="lr",ncol=1)
#ggsave(plot=Fig3, "figs/Enzymesrep/Fig3.png",width=22,height=14)
Fig4<-plot_grid(AllEnzMutsplot,AllEnzCosineplot,AllEnzCosineplotAPO,prop_plot_enzymes_APO,align="v", axis="lr",ncol=1)
#ggsave(plot=Fig4, "figs/Enzymesrep/Fig4.png",width=22,height=19)






AllEnzdiffplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), diff)) +  
  geom_point(alpha = 0.25)+geom_boxplot(alpha=0.5,outlier.shape = NA) + 
  scale_colour_manual(values=c("blue", "red"))+
  #ylim(-5,5)+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  ylab("Difference between Library and WGS TMB estimate")+xlab("Library")
#ggsave(plot=AllEnzdiffplot, "figs/Enzymesrep/AllEnzdiffplot.png",width=22, height =8)


AllEnzdiffplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), diff,col=APOBECHyper)) + geom_boxplot(outlier.shape = NA) + 
  #geom_point(alpha = 0.05,aes(col=APOBECHyper))+
  scale_colour_manual(values=c("blue", "red"))+
  #ylim(-5,5)+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  ylab("Difference between Library and WGS TMB estimate")+xlab("Library")
#ggsave(plot=AllEnzdiffplotAPO, "figs/Enzymesrep/AllEnzdiffplotAPO.png",width=22, height =8)


AllEnzabsdiffplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), abs(diff)))  + 
  geom_jitter(alpha = 0.25,width=0.2,size=5)+geom_boxplot(alpha=0.5,outlier.shape = NA) + ylim(c(0,10))+
  scale_colour_manual(values=c("blue", "red"))+
  #ylim(-5,5)+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  ylab("Absolute value of \ndifference between Library \nand WGS TMB estimate")+xlab("Library")
#ggsave(plot=AllEnzabsdiffplot, "figs/Enzymesrep/AllEnzabsdiffplot.png",width=22, height =8)


AllEnzabsdiffplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), abs(diff),col=APOBECHyper)) + 
  geom_point(alpha = 0.05,position=position_dodge(width=0.75))+geom_boxplot(alpha=0.5,outlier.shape = NA) + 
  scale_colour_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  ylab("Absolute value of \ndifference between Library \nand WGS TMB estimate")+xlab("Library")
#ggsave(plot=AllEnzabsdiffplotAPO, "figs/Enzymesrep/AllEnzabsdiffplotAPO.png",width=22, height =8)

AllEnzArctanplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), arctan)) +
  geom_jitter(alpha = 0.25,width=0.2,size=5)+ geom_boxplot(outlier.shape = NA,alpha=0.5) + 
  scale_colour_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  ylab("Estimated bias of library θ \n(Li,Lo, Nature Scientific \nReports 2021)")+xlab("Library")
#ggsave(plot=AllEnzArctanplot, "figs/Enzymesrep/AllEnzArctanplot.png",width=22, height =8)

AllEnzArctanplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), arctan,col=APOBECHyper)) + geom_boxplot(outlier.shape = NA) + 
  #geom_point(alpha = 0.05,aes(col=APOBECHyper))+
  scale_colour_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  ylab("Estimated bias of library θ \n(Li,Lo, Nature Scientific Reports 2021)")+xlab("Library")
#ggsave(plot=AllEnzArctanplotAPO, "figs/Enzymesrep/AllEnzArctanplotAPO.png",width=22, height =8)
libs<-c("F1","SpeI")
MLcors<-readRDS("MLcors.RDS")



AllEnzInfo<-readRDS("AllEnzInfo.RDS")
AllEnzInfobound<-bind_rows(AllEnzInfo)
ssEnz<-AllEnzInfobound[AllEnzInfobound$Library %in% libs,]
ssEnz
ssEnz$Library<-factor(ssEnz$Library, libs)
corplotLibsF1<-ggplot(ssEnz, aes(mutratewgs,mutrategbs,col=Library))+
  geom_point(alpha=0.25,size=5)+geom_abline(aes(intercept=0, slope=1))+
  xlab("Mutation rate (mutations/Mb)")+ylab("Mutation rate estimate \nfrom Library (mutations/Mb)")+
  theme(axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 25), legend.title = element_text(size = 30),strip.text = element_text(size = 30))+
  guides(colour=guide_legend(override.aes = list(alpha=1)))+facet_wrap(~Library)
corplotLibsF1
corplotLibsF1lim<-ggplot(ssEnz, aes(mutratewgs,mutrategbs,col=Library))+xlim(c(0,10))+ylim(c(0,10))+
  geom_point(alpha=0.25,size=5)+geom_abline(aes(intercept=0, slope=1))+
  xlab("Mutation rate (mutations/Mb)")+ylab("Mutation rate estimate \nfrom Library (mutations/Mb)")+
  theme(axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
        axis.text.x=element_text(vjust=0.5),
        legend.text = element_text(size = 25), legend.title = element_text(size = 30),strip.text = element_text(size = 30))+
  guides(colour=guide_legend(override.aes = list(alpha=1)))+facet_wrap(~Library)
corplotLibsF1lim

MLcors<-readRDS("MLcors.RDS")


MLcors$cor.adj.p<-p.adjust(MLcors$cor.p.val)
MLcors$F.adj.p<-p.adjust(MLcors$Fp.val)
MLcors$pct.genom<-MLcors$gbscov*100/wgscov
MLcors$pct.genom
CorPlot<-ggplot(MLcors, aes(pct.genom,cor)) + geom_point(size=5)+
  geom_point(data=MLcors[MLcors$Library=="F1",], aes(pct.genom,cor),col="red", size=5)+
  geom_point(data=MLcors[MLcors$Library=="SpeI",], aes(pct.genom,cor),col="blue", size=5)+
  theme(axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
        axis.text.x=element_text(vjust=0.5),
        legend.text = element_text(size = 25), legend.title = element_text(size = 30),strip.text = element_text(size = 30))+
  ylab("Correlation to WGS \nTMB estimate")+
  xlab("Percentage genome captured \n(log2 scale)")+ylim(c(0,1))+scale_x_continuous(trans = 'log2')
CorPlot
#MLcors<-MLcors[MLcors$cor>0.9,]
#ggsave(plot=CorPlot, "figs/Enzymesrep/CorPlot.png",width=30,height=14)
arctan_mod<-lm(log2(median_arctan) ~ log2(pct.genom),MLcors)
eq_arctan<-summary(arctan_mod)
AdjR2_arctan<-eq_arctan$adj.r.squared
r_arctan<-eq_arctan$coefficients[2,1]
p_arctan<-eq_arctan$coefficients[2,4]
annot<-paste0("Adjusted R^2=",signif(AdjR2_arctan,3),
              "\np=",signif(p_arctan,3))
MLcors$estim_arctan<-log2(MLcors$median_arctan)-eq_arctan$coefficients[2,1]*log2(MLcors$gbscov)-eq_arctan$coefficients[1,1]
d<-seq(from=min(MLcors$pct.genom),to=max(MLcors$pct.genom),length.out = 10000)
df<-data.frame(pct.genom=d)
dfnew<-data.frame(pct.genom=d,2^predict(arctan_mod,newdata=df,interval='confidence'))
Arctanplot<-ggplot(MLcors,aes(pct.genom,median_arctan))+  
  geom_line(data=dfnew,(aes(pct.genom, fit)),col="blue")+
  geom_line(data=dfnew,(aes(pct.genom, lwr)),col="blue",linetype="dotted")+
  geom_line(data=dfnew,(aes(pct.genom, upr)),col="blue",linetype="dotted")+
  geom_point(size=3)+
  geom_point(data=MLcors[MLcors$Library=="F1",], aes(pct.genom,median_arctan),col="red", size=5)+
  geom_text(aes(label=ifelse(Library=="F1",as.character(Library),'')),size=10,hjust=-0.3,col="red")+
  geom_point(data=MLcors[MLcors$Library=="SpeI",], aes(pct.genom,median_arctan),col="blue", size=5)+
  geom_text(aes(label=ifelse(Library=="SpeI",as.character(Library),'')),size=10,hjust=-0.3,col="blue")+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  geom_text(aes(label=annot,x=25,y=0.25),size=5)+
  ylab("Median estimated \nbias of library θ")+xlab("Library Coverage (log2 scale)")+
  scale_x_continuous(trans='log2')
Arctanplot
#ggsave(plot=Arctanplot, "figs/Enzymes/Arctanplot.png",width =10,height=10)

diff_mod<-lm(log2(median_diff) ~ log2(pct.genom),MLcors)
eq_diff<-summary(diff_mod)
AdjR2_diff<-eq_diff$adj.r.squared
r_diff<-eq_diff$coefficients[2,1]
p_diff<-eq_diff$coefficients[2,4]
annot_diff<-paste0("Adjusted R^2=",signif(AdjR2_diff,3),
              "\np=",signif(p_diff,3))
MLcors$estim_diff<-log2(MLcors$median_diff)-eq_diff$coefficients[2,1]*log2(MLcors$pct.genom)-eq_diff$coefficients[1,1]
d<-seq(from=min(MLcors$pct.genom),to=max(MLcors$pct.genom),length.out = 10000)
df<-data.frame(pct.genom=d)
dfnew<-data.frame(pct.genom=d,2^predict(diff_mod,newdata=df,interval='confidence'))
diffplot<-ggplot(MLcors,aes(pct.genom,median_diff))+
  geom_line(data=dfnew,(aes(pct.genom, fit)),col="blue")+
  geom_line(data=dfnew,(aes(pct.genom, lwr)),col="blue",linetype="dotted")+
  geom_line(data=dfnew,(aes(pct.genom, upr)),col="blue",linetype="dotted")+
  geom_point(size=5)+
  geom_point(data=MLcors[MLcors$Library=="F1",], aes(pct.genom,median_diff),col="red",size=8)+
  geom_point(data=MLcors[MLcors$Library=="SpeI",], aes(pct.genom,median_diff),col="blue",size=8)+
  geom_text(aes(label=ifelse(Library=="F1",as.character(Library),'')),size=10,hjust=-0.3,col="red")+
  geom_text(aes(label=ifelse(Library=="SpeI",as.character(Library),'')),size=10,hjust=-0.3,col="blue")+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  geom_text(aes(label=annot_diff,x=25,y=0.6),size=10)+
  ylab("Median absolute difference \nbetween WGS TMB and \n Library TMB")+xlab("Percentage genome captured \n(log2 scale)")+ylim(c(0,1))+
  scale_x_continuous(trans = 'log2')
diffplot

plot_grid(CorPlot,corplotLibsF1,corplotLibsF1lim,diffplot,ncol=1, align="v",axis="lr")
#ggsave("figs/Enzymesrep/Fig5.png",width=18,height=24)
?plot_grid
gp<-plot_grid(AllEnzCosineplot+theme(axis.title = element_text(size=50), axis.text = element_text(size=45),
                                 axis.text.x = element_text(angle=45,hjust=1,vjust=1)),
          CorPlot+theme(axis.title = element_text(size=50), axis.text = element_text(size=45)),
          corplotLibsF1lim+theme(axis.title = element_text(size=50), axis.text = element_text(size=45),
                                 legend.text = element_text(size = 40), 
                                 legend.title = element_text(size = 45),strip.text = element_text(size = 40)),
          diffplot+theme(axis.title = element_text(size=50), axis.text = element_text(size=45)),ncol=1, align="v",axis="lr",scale=0.9)
ggsave(plot=gp,"figs/Enzymesrep/posterplot.png",width=25,height=40)
diffplotnolog<-ggplot(MLcors,aes(pct.genom,median_diff))+
  geom_line(data=dfnew,(aes(pct.genom, fit)),col="blue")+
  geom_line(data=dfnew,(aes(pct.genom, lwr)),col="blue",linetype="dotted")+
  geom_line(data=dfnew,(aes(pct.genom, upr)),col="blue",linetype="dotted")+
   geom_point(size=1)+
  # geom_point(data=MLcors[MLcors$Library=="F1",], aes(pct.genom,median_diff),col="red",size=4)+
  # geom_point(data=MLcors[MLcors$Library=="SpeI",], aes(pct.genom,median_diff),col="blue",size=4)+
  # geom_text(aes(label=ifelse(Library=="F1",as.character(Library),'')),size=10,hjust=-0.3,col="red")+
  # geom_text(aes(label=ifelse(Library=="SpeI",as.character(Library),'')),size=10,hjust=-0.3,col="blue")+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  geom_text(aes(label=annot_diff,x=25,y=0.6),size=10)+
  ylab("Median absolute difference \nbetween WGS TMB and \n Library TMB")+xlab("Percentage genome captured")+ylim(c(0,1))
diffplotnolog
error<-plot_grid(Arctanplot,diffplot,nrow = 1)

plot_grid(corplotLibs,error,ncol=1)
AllEnzInfo<-readRDS("AllEnzInfofilt.RDS")
AllEnzInfobound<-bind_rows(AllEnzInfo)
libtest<-MLcors$Library[log2(MLcors$gbscov)<21 & log2(MLcors$gbscov)>20.5]
libtest<-AllEnzInfobound[AllEnzInfobound$Library %in% libtest,]

libtestdiff<-ggplot(libtest, aes(reorder(Library,nmutsgbs), abs(diff)))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.1,alpha=0.1)+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  xlab("Library")+
  ylab("Absolute difference \n between WGS TMB and \n Library TMB")

libtestarctan<-ggplot(libtest, aes(reorder(Library,nmutsgbs), arctan))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.1,alpha=0.1)+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 30), axis.text =  element_text(size = 25),
    axis.text.x=element_text(vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  xlab("Library")+
  ylab("Estimated bias of library θ \n(Li,Lo, Nature Scientific \nReports 2021, log2 scale)")
libtests<-plot_grid(diffplot,Arctanplot,libtestdiff,libtestarctan,nrow = 2,ncol=2,align="v")
libtests
#ggsave("figs/Enzymesrep/Fig6.png",width=18,height=24)
CorPlotgrid<-CorPlot+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
AllEnzabsdiffplotgrid<-AllEnzabsdiffplot+theme(axis.text.x = element_blank(),axis.title.x = element_blank())+ylim(c(0,10))
AllEnzMutsplotgrid<-AllEnzMutsplot+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
Fig6<-plot_grid(CorPlotgrid, AllEnzabsdiffplot, align="v",axis="lr",ncol=1)
#ggsave(plot=Fig6, "figs/Enzymesrep/Fig6.png",width=18,height=25)
sublist<-subset(MLcors, pct.genom < 1.3 & median_diff <0.125)
