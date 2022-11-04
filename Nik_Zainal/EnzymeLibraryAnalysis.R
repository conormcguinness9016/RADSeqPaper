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
library(parallel)
library(ggrepel)
setwd("Nik_Zainal/")

#REfraglibs<-readRDS("REfraglibs.RDS")
#REfraglibs<-REfraglibs[lapply(REfraglibs, function(x){nrow(x)}) >0]
AllEnzInfo<-readRDS("AllEnzInfo.RDS")
dfprops<-readRDS("dfprops.RDS")
dfprops

print("read in enzyme information")





library(GenomicFeatures)
f1<-scan("FoundationOneGenes.txt", what = "", sep = " ")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exons<-exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
library(biomaRt)
ensembl<-useEnsembl(biomart="ensembl",GRCh=37, dataset = "hsapiens_gene_ensembl")
library(stringr)
geneinfo<-getBM(ensembl,filters = "hgnc_symbol", attributes = c("hgnc_id", "hgnc_symbol", "entrezgene_id"), values = f1)
geneinfo$entrezgene_id <- as.character(geneinfo$entrezgene_id)
f1ranges<-exons[which(names(exons) %in% geneinfo$entrezgene_id)]
covs<-mclapply(seq(REfraglibs),function(x){
   if(nrow(REfraglibs[[x]])>0){
   print(x)
   RE.gr<-makeGRangesFromDataFrame(REfraglibs[[x]])
   gbscov<-sum(width(GenomicRanges::reduce(RE.gr)))
   data.frame(Library=REfraglibs[[x]]$Library[1],gbscov)}
 }, mc.cores = 10)

 covs<-bind_rows(covs)
 covs<-bind_rows(covs,data.frame(Library="F1", gbscov=sum(sum(width(GenomicRanges::reduce(f1ranges))))))
 covs$Library<-str_extract(covs$Library, "\\w{1,}\\p{Uppercase}")
 covs$Library[length(covs$Library)]<-"F1"
 covs<-covs[!duplicated(covs$Library),]
 covs<-covs[order(covs$gbscov),]
 covs
saveRDS(covs,"covs.RDS")
covs<-readRDS("covs.RDS")

enz<-names(AllEnzInfo)
for(i in seq(AllEnzInfo)){
  AllEnzInfo[[i]]$APOBECHyper<-dfprops$APOBECHyper[match(AllEnzInfo[[i]]$Sample, dfprops$Var2)]
  AllEnzInfo[[i]]$diff<-AllEnzInfo[[i]]$mutrategbs-AllEnzInfo[[i]]$mutratewgs
  AllEnzInfo[[i]]$mutclasswgs<-NA
  AllEnzInfo[[i]]$mutclassgbs<-NA
  AllEnzInfo[[i]]$Library<-rep(enz[i],nrow(AllEnzInfo[[i]]))
  AllEnzInfo[[i]]$mutclasswgs[AllEnzInfo[[i]]$mutratewgs>10]<-"High (>10muts/Mb)"
  AllEnzInfo[[i]]$mutclasswgs[AllEnzInfo[[i]]$mutratewgs<10 & AllEnzInfo[[i]]$mutratewgs>5]<-"Medium (5-10muts/Mb)"
  AllEnzInfo[[i]]$mutclasswgs[AllEnzInfo[[i]]$mutratewgs<5]<-"Low (<5muts/Mb)"
  AllEnzInfo[[i]]$mutclassgbs[AllEnzInfo[[i]]$mutrategbs>10]<-"High (>10muts/Mb)"
  AllEnzInfo[[i]]$mutclassgbs[AllEnzInfo[[i]]$mutrategbs<10 & AllEnzInfo[[i]]$mutrategbs>5]<-"Medium (5-10muts/Mb)"
  AllEnzInfo[[i]]$mutclassgbs[AllEnzInfo[[i]]$mutrategbs<5]<-"Low (<5muts/Mb)"
  AllEnzInfo[[i]]$arctan<-abs((pi/4)-atan2(AllEnzInfo[[i]]$mutratewgs,AllEnzInfo[[i]]$mutrategbs))
}
AllEnzInfo
saveRDS(AllEnzInfo, "AllEnzInfo.RDS")
F1cov<-sum(sum(width(GenomicRanges::reduce(f1ranges))))
wgscov<-sum(seqlengths(Hsapiens)[1:24])

MLcors<-lapply(seq(length(AllEnzInfo)), function(x){
  print(x)
  df<-AllEnzInfo[[x]]
  Library<-df$Library[1]
  gbscov<-covs$gbscov[covs$Library %in% Library]
  mean_nmutsgbs<-mean(df$nmutsgbs)
  if(mean_nmutsgbs>2){
    median_nmutsgbs<-median(df$nmutsgbs)
    median_diff<-median(abs(df$diff))
    median_arctan<-median(df$arctan)
    mean_arctan<-mean(df$arctan)
    mean_diff<-mean(abs(df$diff))
    cor<-cor.test(df$mutrategbs, df$mutratewgs,method="pearson")
    var<-var.test(df$mutrategbs, df$mutratewgs)
    t<-t.test(df$mutrategbs,df$mutratewgs,paired=TRUE)
    data.frame(Mean_mutations_captured=mean_nmutsgbs,median_nmutsgbs,
               median_diff,mean_diff,median_arctan,mean_arctan,
               gbscov,
               cor=cor$estimate,corlwr=cor$conf.int[1],
               corupr=cor$conf.int[2],
               cor.p.val=cor$p.value,
               Library=df$Library[1],
               Fvar=var$statistic,
               Fp.val=var$p.value,
               Flwr=var$conf.int[1],
               Fupr=var$conf.int[2],
               t=t$statistic,
               t.pval=t$p.value)
 }
  })
MLcors<-bind_rows(MLcors)
MLcors$cor.adj.p<-p.adjust(MLcors$cor.p.val)
MLcors$F.adj.p<-p.adjust(MLcors$Fp.val)
saveRDS(MLcors,"MLcors.RDS")

CorPlot<-ggplot(MLcors, aes(log2(gbscov),cor)) + geom_point(size=3)+
  geom_point(data=MLcors[MLcors$Library=="F1",], aes(log2(gbscov),cor),col="red", size=5)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5), axis.title = element_text(size=15),
        axis.text = element_text(size=10))+ylab("Correlation to whole \ngenome TMB (mutations/Mb) estimate")+
  xlab("Library Coverage(log2 scale)")+ylim(c(0,1))
CorPlot
ggsave(plot=CorPlot, "figs/Enzymesrep/CorPlot.png",width=22,height=14)
#MLcors<-MLcors[MLcors$cor>0.9,]
AllEnzInfo<-AllEnzInfo[names(AllEnzInfo) %in% MLcors$Library]
AllEnzInfo
saveRDS(AllEnzInfo, "AllEnzInfofilt.RDS")
accept<-0.95
takeoff<-lapply(seq(AllEnzInfo), function(x){
  df<-AllEnzInfo[[x]]
  Ah_tests<-table(df$APOBECHyper)[2]
  At_tests<-table(df$APOBECHyper)[1]
  All_tests<-sum(table(df$APOBECHyper))
  test<-table(df$cosimdiag>0.9, df$APOBECHyper)
  TF<-data.frame(test)
  TF<-subset(data.frame(test), Var1==TRUE)
  if( nrow(TF)==0){
    TF<-data.frame(Var1=c(TRUE,TRUE), Var2 = c("APOBEC_hyper", "APOBEC_typical"), Freq=c(0,0))
  }
  Ah_success<-TF$Freq[TF$Var2=="APOBEC_hyper"]
  At_success<-TF$Freq[TF$Var2=="APOBEC_typical"]
  All_success<-sum(TF$Freq)
  Ah_prop<-Ah_success/Ah_tests
  At_prop<-At_success/At_tests
  All_prop<-All_success/All_tests
  Ah_success
  binomAh<-binom.test(Ah_success, Ah_tests, accept, alternative = "greater")
  resAh<-data.frame(APOBEChyper="APOBEC_hyper",
                    Library=names(AllEnzInfo[x]),
                    mean_nmutsgbs=mean(AllEnzInfo[[x]]$nmutsgbs),
                    p.val=binomAh$p.value,
                    lwr_95=binomAh$conf.int[2],
                    upr_95=binomAh$conf.int[1],
                    estimate=as.numeric(binomAh$estimate),
                    tests=as.numeric(binomAh$parameter),
                    prop=as.numeric(Ah_prop))
  binomAt<-binom.test(x=At_success, n=At_tests,accept, alternative = "greater")
  resAt<-data.frame(APOBEChyper="APOBEC_typical",
                    Library=names(AllEnzInfo[x]),
                    mean_nmutsgbs=mean(AllEnzInfo[[x]]$nmutsgbs),
                    p.val=binomAt$p.value,
                    lwr_95=binomAt$conf.int[1],
                    upr_95=binomAt$conf.int[2],
                    estimate=as.numeric(binomAt$estimate),
                    tests=as.numeric(binomAt$parameter),
                    prop=as.numeric(At_prop))
  binomAll<-binom.test(x=All_success, n=All_tests,accept, alternative = "greater")
  resAll<-data.frame(APOBEChyper="All_samples",
                     Library=names(AllEnzInfo[x]),
                     mean_nmutsgbs=mean(AllEnzInfo[[x]]$nmutsgbs),
                     p.val=binomAll$p.value,
                     lwr_95=binomAll$conf.int[1],
                     upr_95=binomAll$conf.int[2],
                     estimate=as.numeric(binomAll$estimate),
                     tests=as.numeric(binomAll$parameter),
                     prop=as.numeric(All_prop))
  bind_rows(resAh, resAt,resAll)
}
)

takeoff<-bind_rows(takeoff)
takeoff$adj.pval<-p.adjust(takeoff$p.val)
takeoff$sig<-""
takeoff$sig[takeoff$adj.pval<0.05]<-"*"
takeoff$APOBEChyper<-fct_relevel(takeoff$APOBEChyper, c("APOBEC_typical","APOBEC_hyper", "All_samples"))
saveRDS(takeoff,"takeoff.RDS")
takeoff_All<-subset(takeoff, APOBEChyper == "All_samples")
prop_plot_enzymes_all<-ggplot(takeoff_All, aes(reorder(Library, mean_nmutsgbs), prop)) + geom_bar(stat="identity", position="dodge")+#geom_point()+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  #geom_bar(stat="identity", position ="dodge")+
  # geom_text(aes(x=first_sig_Ah,y=1.02, colour =APOBEChyper[1], label="*"), size =10, show.legend = FALSE)+
  geom_text(aes(x=reorder(Library, mean_nmutsgbs),y=1.02, label=sig), size =5,show.legend = FALSE)+
  xlab("Library")+
  ylab("Proportion of samples \nresulting in >90% \ncosine similarity") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 10),
    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"))
prop_plot_enzymes_all
dir.create("figs/Enzymes/")
ggsave(plot=prop_plot_enzymes_all,"figs/Enzymes/Proportion_plot_All_Enzymes.png",width=20,height=8)

takeoff_APO<-subset(takeoff, APOBEChyper != "All_samples")
prop_plot_enzymes_APO<-ggplot(takeoff_APO, aes(reorder(Library, mean_nmutsgbs), prop,fill=APOBEChyper)) + geom_bar(stat="identity", position="dodge")+#geom_point()+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  #geom_bar(stat="identity", position ="dodge")+
  # geom_text(aes(x=first_sig_Ah,y=1.02, colour =APOBEChyper[1], label="*"), size =10, show.legend = FALSE)+
  geom_text(aes(x=reorder(Library, mean_nmutsgbs),y=1.02, label=sig, col =APOBEChyper), size =10,show.legend = FALSE)+
  xlab("Library")+
  ylab("Proportion of samples \nresulting in >90% \ncosine similarity") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 10),
    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"))
prop_plot_enzymes_APO
ggsave(plot=prop_plot_enzymes_APO,"figs/Enzymes/Proportion_plot_Apo_Enzymes.png",width=22,height=8)


AllEnzInfobound<-bind_rows(AllEnzInfo)
AllEnzCosineplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), cosimdiag)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.025)+
  scale_colour_manual(values=c("blue", "red"))+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5), axis.title = element_text(size=15),
        axis.text = element_text(size=10))+ylab("Cosine similarity \n to WGS profile")+xlab("Library")
ggsave(plot=AllEnzCosineplot, "figs/Enzymes/AllEnzCosineplot.png",width=22,height=8)
AllEnzCosineplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), cosimdiag,col=APOBECHyper)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75),alpha = 0.05)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 10),
    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"),
         override.aes=list(alpha=1))+ylab("Cosine similarity \n to WGS profile")+xlab("Library")
AllEnzMutsplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), nmutsgbs)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.025)+
  scale_colour_manual(values=c("blue", "red"))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),axis.text=element_text( size =10), axis.title = element_text(size=15))+
  ylab("Number of mutations \ndetected (log10 scale))")+xlab("Library")+
  scale_y_continuous(trans='log10')
ggsave(plot=AllEnzMutsplot, "figs/Enzymes/AllEnzMutsplot.png",width=22, height =8)
ggsave(plot=AllEnzCosineplotAPO, "figs/Enzymes/AllEnzCosineplotAPO.png",width=22,height=8)
Supp_fig1<-plot_grid(AllEnzMutsplot,AllEnzCosineplot,AllEnzCosineplotAPO,prop_plot_enzymes_APO,align="v", axis="lr",ncol=1)
ggsave(plot=Supp_fig1, "figs/Enzymes/Supp_fig1.png",width=22,height=18)




AllEnzdiffplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), diff)) +  
  geom_point(alpha = 0.05)+geom_boxplot(alpha=0.5,outlier.shape = NA) + 
  scale_colour_manual(values=c("blue", "red"))+
  #ylim(-5,5)+
  theme(axis.text=element_text( size =10), axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),axis.title = element_text(size=15))+
  ylab("Difference between Library and WGS TMB estimate")+xlab("Library")
ggsave(plot=AllEnzdiffplot, "figs/Enzymes/AllEnzdiffplot.png",width=22, height =8)


AllEnzdiffplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), diff,col=APOBECHyper)) + geom_boxplot(outlier.shape = NA) + 
  #geom_point(alpha = 0.05,aes(col=APOBECHyper))+
  scale_colour_manual(values=c("blue", "red"))+
  #ylim(-5,5)+
  theme(axis.text=element_text( size =10),axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5), axis.title = element_text(size=15))+
  ylab("Difference between Library and WGS TMB estimate")+xlab("Library")
ggsave(plot=AllEnzdiffplotAPO, "figs/Enzymes/AllEnzdiffplotAPO.png",width=22, height =8)


AllEnzabsdiffplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), abs(diff)))  + 
  geom_point(alpha = 0.05)+geom_boxplot(alpha=0.5,outlier.shape = NA) + 
  scale_colour_manual(values=c("blue", "red"))+
  #ylim(-5,5)+
  theme(axis.text=element_text( size =10),axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5), axis.title = element_text(size=15))+
  ylab("Absolute value of \ndifference between Library and WGS TMB estimate")+xlab("Library")
ggsave(plot=AllEnzabsdiffplot, "figs/Enzymes/AllEnzabsdiffplot.png",width=22, height =8)


AllEnzabsdiffplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), abs(diff),col=APOBECHyper)) + 
  geom_point(alpha = 0.05,position=position_dodge(width=0.75))+geom_boxplot(alpha=0.5,outlier.shape = NA) + 
  scale_colour_manual(values=c("blue", "red"))+
  theme(axis.text=element_text(size =10),axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5), axis.title = element_text(size=15))+
  ylab("Absolute value of \ndifference between Library and WGS TMB estimate")+xlab("Library")
ggsave(plot=AllEnzabsdiffplotAPO, "figs/Enzymes/AllEnzabsdiffplotAPO.png",width=22, height =8)

AllEnzArctanplot<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), arctan)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.05)+
  scale_colour_manual(values=c("blue", "red"))+
  theme(axis.text=element_text( size =10), axis.title = element_text(size=15),axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+
  ylab("Estimated bias of library θ \n(Li,Lo, Nature Scientific Reports 2021)")+xlab("Library")
ggsave(plot=AllEnzArctanplot, "figs/Enzymes/AllEnzArctanplot.png",width=22, height =8)

AllEnzArctanplotAPO<-ggplot(AllEnzInfobound, aes(reorder(Library, nmutsgbs), arctan,col=APOBECHyper)) + geom_boxplot(outlier.shape = NA) + 
  #geom_point(alpha = 0.05,aes(col=APOBECHyper))+
  scale_colour_manual(values=c("blue", "red"))+
  theme(axis.text=element_text( size =10), axis.title = element_text(size=15),axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+
  ylab("Estimated bias of library θ \n(Li,Lo, Nature Scientific Reports 2021)")+xlab("Library")
ggsave(plot=AllEnzArctanplotAPO, "figs/Enzymes/AllEnzArctanplotAPO.png",width=22, height =8)
CorPlotgrid<-CorPlot+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
AllEnzabsdiffplotgrid<-AllEnzabsdiffplot+theme(axis.text.x = element_blank(),axis.title.x = element_blank())+ylim(c(0,10))
AllEnzMutsplotgrid<-AllEnzMutsplot+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
Supp_fig_3<-plot_grid(AllEnzArctanplot, AllEnzabsdiffplot+ylim(c(0,10)), align="v",axis="lr",ncol=1)
ggsave(plot=Supp_fig_3, "figs/Enzymes/Supp_fig3.png",width=22,height=12)



