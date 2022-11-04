setwd("/Volumes/archive/cancergeneticslab/ConorM/RADSeqPaper/")

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
library(data.table)
library(parallel)
library(ggsignif)
library(gridExtra)
library(gtable)
library(grid)

print("libraries loaded")
#frags<-readRDS("REfraglibs.rds")
print("loaded frags")

print("")
setwd("Nik_Zainal")
BRCAgenoms<-read.delim("geninfo560genomes.txt")
BRCASigProf<-data.frame(Project=rep("SigProf",nrow(BRCAgenoms)),
                        Sample=BRCAgenoms$Sample,
                        ID=rep(".",nrow(BRCAgenoms)),
                        Genome=rep("GRCh37",nrow(BRCAgenoms)),
                        mut_type=rep("SNP",nrow(BRCAgenoms)),
                        chrom=BRCAgenoms$Chrom,
                        pos_start=BRCAgenoms$Pos,
                        pos_end=BRCAgenoms$Pos,
                        ref=BRCAgenoms$Ref,
                        alt=BRCAgenoms$Alt,
                        Type=rep("SOMATIC",nrow(BRCAgenoms)))
dir.create("SigProf")
write.table(BRCASigProf, file="SigProf/BRCANZ.txt",sep=" ",row.names = FALSE,quote = FALSE,)

BRCAgenoms$start<-BRCAgenoms$Pos
BRCAgenoms$end<-BRCAgenoms$Pos
BRCAgenoms$Chrom<-unlist(lapply(BRCAgenoms$Chrom, function(x){
  paste0("chr", x)
}))
BRCAgenoms.gr<-makeGRangesFromDataFrame(BRCAgenoms, keep.extra.columns = TRUE)
BRCA_vr = VRanges(
  seqnames = seqnames(BRCAgenoms.gr),
  ranges = ranges(BRCAgenoms.gr),
  ref = BRCAgenoms.gr$Ref,
  alt = BRCAgenoms.gr$Alt,
  sampleNames = BRCAgenoms.gr$Sample,
  study = "overlap")
BRCA_motifs<-mutationContext(BRCA_vr, BSgenome.Hsapiens.UCSC.hg19)
saveRDS(BRCA_motifs, "BRCA_motifs.RDS")

BRCA_motifs<-readRDS("BRCA_motifs.RDS")

btab<-as.data.table(BRCA_motifs)

btab$Context<-"Other"
btab$Context[btab$context =="T.A" & btab$alteration=="CT"]<-"APOBEC"
btab$Context[btab$context =="T.T" & btab$alteration=="CT"]<-"APOBEC"
btab$Context[btab$context =="T.A" & btab$alteration=="CG"]<-"APOBEC"
btab$Context[btab$context =="T.T" & btab$alteration=="CG"]<-"APOBEC"
BRCA_mm = motifMatrix(BRCA_motifs, group = "sampleNames", normalize = TRUE)
BRCAmatrix = motifMatrix(BRCA_motifs, group = "sampleNames", normalize = FALSE)
wgscov<-sum(seqlengths(Hsapiens)[1:24])
nmutswgs<-colSums(BRCAmatrix)

dfprops<-data.frame(prop.table(table(btab$Context, btab$sampleNames),2))
dfprops<-subset(dfprops, Var1 == "APOBEC")
dfprops$nmuts<-nmutswgs
dfprops$TMB[dfprops$nmuts>30000]<-"High (>10muts/Mb)"
dfprops$TMB[dfprops$nmuts<30000 & dfprops$nmuts>15000]<-"Medium (5-10muts/Mb)"
dfprops$TMB[dfprops$nmuts<15000]<-"Low (<5muts/Mb)"
dfprops$APOBECHyper[dfprops$Freq>0.4]<-"APOBEC_hyper"
dfprops$APOBECHyper[dfprops$Freq<0.4]<-"APOBEC_typical"

dfprops$APOBECHyper<-fct_relevel(dfprops$APOBECHyper,c("APOBEC_typical", "APOBEC_hyper"))
dfprops$TMB<-as.factor(dfprops$TMB)
dfprops$TMB<-fct_relevel(dfprops$TMB, c("Low (<5muts/Mb)", "Medium (5-10muts/Mb)", "High (>10muts/Mb)"))
saveRDS(dfprops, "dfprops.RDS")
dfprops<-readRDS("dfprops.RDS")
fish<-fisher.test(table(dfprops$TMB, dfprops$APOBECHyper))

dir.create("figs")
ggplot(dfprops, aes(TMB, fill = APOBECHyper))+ geom_bar()+
  geom_text(aes(label=paste0("Fisher's exact \np=",signif(fish$p.value, 3)), y=Inf, x=3), size =5)+
  ylab("Count")+xlab("Tumour mutation burden")+
  coord_cartesian(clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill = guide_legend(title="APOBEC \nsignature status"))

ggsave("figs/TMB_by_APOBEC.png",width=10,height=5)

ggplot(dfprops, aes(TMB, fill = APOBECHyper))+ geom_bar(position = "fill")+
  geom_text(aes(label=paste0("Fisher's exact \np=",signif(fish$p.value, 3)), y=Inf, x=3), size =5)+
  ylab("Proportion")+xlab("Tumour mutation burden")+
  coord_cartesian(clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill = guide_legend(title="APOBEC \nsignature status"))
ggsave("figs/TMB_by_APOBEC_proportion.png",width=10,height=5)



ggplot(dfprops, aes(APOBECHyper, fill = TMB))+ geom_bar()+
  geom_text(aes(label=paste0("Fisher's exact \np=",signif(fish$p.value, 3)), y=Inf, x=2), size =5)+
  ylab("Count")+xlab("APOBEC status")+
  coord_cartesian(clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  scale_fill_manual(values=c("pink","plum", "purple"))+
  guides(fill=guide_legend(title="Tumour mutation \nburden"))

ggsave("figs/APOBEC_by_TMB.png",width=10,height=5)
ggplot(dfprops, aes(APOBECHyper, fill = TMB))+ geom_bar(position="fill")+
  geom_text(aes(label=paste0("Fisher's exact \np=",signif(fish$p.value, 3)), y=Inf, x=2), size =5)+
  ylab("Proportion")+xlab("APOBEC status")+
  coord_cartesian(clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  scale_fill_manual(values=c("pink","plum", "purple"))+
  guides(fill=guide_legend(title="Tumour mutation \nburden"))

ggsave("figs/APOBEC_by_TMB_proportion.png",width=10,height=5)
w<-wilcox.test(nmuts ~ APOBECHyper, dfprops)
w
##medians
aggregate(.~APOBECHyper, dfprops, median)
##sample sizes
table(dfprops$APOBECHyper)
TMBdf<-data.frame(xmin=rep(-Inf,3),
                  xmax=rep(Inf,3),
                  ymin=c(-Inf,15000,30000),
                  ymax=c(15000,30000,Inf),
                  TMB=c("Low (<5muts/Mb)","Medium (5-10muts/Mb)","High (>10muts/Mb)")
)
tmbtab<-table(TMB=dfprops$TMB, APOBEC=dfprops$APOBECHyper)
tmbtab<-tmbtab[c(3,2,1),]
colnames(tmbtab)<-c("APOBEC\ntypical", "APOBEC\nhyper")
tab<-table(dfprops$APOBECHyper,dfprops$TMB=="High (>10muts/Mb)")
fisher.test(tab)

tmbtabtg <- tableGrob(tmbtab)
TMB_APO<-ggplot(dfprops, mapping=aes(APOBECHyper, nmuts))+
  geom_rect(data=TMBdf, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=TMB),alpha=0.5,inherit.aes=FALSE)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1)+
  geom_signif(comparisons = list(c("APOBEC_hyper", "APOBEC_typical"), annotation = w$p.value),
              map_signif_level =FALSE, size =1, textsize =5)+
  ylab("Mutation Count (WGS)")+xlab("APOBEC status")+
  coord_cartesian(xlim=c(1,2),ylim=c(0,100000),clip = 'off')+
  theme(#plot.margin = unit(c(1.5,2.5,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.position="right",
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  scale_y_continuous(breaks=c(0,15000,30000,100000))+
  scale_fill_manual(breaks=c("High (>10muts/Mb)","Medium (5-10muts/Mb)","Low (<5muts/Mb)"),values=c("purple","plum", "pink"))+
  guides(fill=guide_legend(title="Tumour mutation \nburden (TMB)"))+
  #annotation_custom(tmbtabtg,xmin=2.18,xmax=4,ymin=45000,ymax=45000)+
  geom_text(aes(label=paste0("Fisher's exact \np=",signif(fish$p.value, 3.5)), y=15000, x=2.85), size =5)
ggsave(plot=TMB_APO,"figs/TMB_byAPOBEC_points.png",width=10,height=5)
TMB_APO
setDT(BRCAgenoms)

##*5*100


randsample<-function(x){
  btab[which(btab$sampleNames %in% names(which(table(btab$sampleNames)>x))),
       .SD[sample(.N, x, replace = FALSE)],
       by = sampleNames]
}

mutinfo<-function(x, mutrate=FALSE){
  overlap_motifs<-makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  overlap_mm = motifMatrix(overlap_motifs, group = "sampleNames", normalize = TRUE)
  overlapmatrix = motifMatrix(overlap_motifs, group = "sampleNames", normalize = FALSE)
  cosmat<-cos_sim_matrix(overlap_mm,BRCA_mm)
  cosimdiag<-diag(cosmat)
  nmutsgbs<-colSums(overlapmatrix)
  df<-data.table(nmutswgs, nmutsgbs, cosimdiag, Sample = names(nmutswgs))
}



now<-Sys.time()
expo<-list()
rando<-NULL
rando<-data.table(mclapply(seq(10,800, by =10), randsample, mc.cores=4))
oss<-mclapply(rando$V1, mutinfo,mc.cores =8)
allrands<-bind_rows(oss)
allrands$APOBECHyper<-dfprops$APOBECHyper[match(allrands$Sample, dfprops$Var2)]

allrands$APOBECrate<-dfprops$Freq[match(allrands$Sample, dfprops$Var2)]
allrands
saveRDS(allrands,"allrands.RDS",)
tests<-lapply(seq(10,800, by =10),function(x){
  d<-subset(allrands, nmutsgbs ==x)
  t<-t.test(cosimdiag ~ APOBECHyper, d)
  data.frame(p.val=t$p.value, nmutsgbs = x)
})
tests<-bind_rows(tests)
tests$p.adj<-p.adjust(tests$p.val, "holm")
ggplot(tests, aes(nmutsgbs, -log10(p.adj))) + geom_line()

dir.create("figs/")
#png"figs/RandomMutsVsCosineSim800.png")
Randplot<-ggplot(allrands, aes(nmutsgbs, cosimdiag))+ geom_point(alpha =0.05)+
  ylab("Cosine similarity \n to WGS profile")+
  xlab("Number of random \n mutations sampled")+
  geom_hline(aes(yintercept=0.90))+
  theme(axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  scale_y_continuous(breaks = c(.25,.5,.75,.9,1))
Randplot
ggsave(plot=Randplot,"figs/RandomMutsVsCosineSim800.png",width=10,height=5)

#dev.off()
scl<-function(x){(x-min(x))/(max(x)-min(x))}
allrands$nmutswgs_sc<-scl(allrands$nmutswgs)

#png"figs/RandomMutsVsCosineSimMLdepicted800.png")
MLplot<-ggplot(allrands, aes(nmutsgbs, cosimdiag, colour =nmutswgs))+
  geom_point(aes(group = Sample), alpha =allrands$nmutswgs_sc)+
  scale_colour_gradient(low = "pink",high="purple", name="Number of \nmutations (WGS)")+
  ylab("Cosine similarity \n to WGS profile")+
  xlab("Number of random \n mutations sampled")+
  theme(axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave(plot=MLplot,"figs/RandomMutsVsCosineSim800MLdepicted800.png",width=10,height=5)
#dev.off()
MLplot
#png"figs/RandomMutsVsCosineSimAPOBEC_rate.png")
Aporate<-ggplot(allrands, aes(nmutsgbs, cosimdiag))+
  geom_point(aes(colour=APOBECrate), alpha = 0.1)+
  ylab("Cosine similarity \n to WGS profile")+
  xlab("Number of random \n mutations sampled")+
  scale_colour_gradient(low="blue", high="red",name="Proportion of \nAPOBEC mutations")+
  theme(axis.title = element_text(size =18), axis.text =  element_text(size = 15),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave(plot=Aporate,"figs/RandomMutsVsCosineSimAPOBEC_rate.png",width=10,height=5)
#dev.off()
Aporate
#png"figs/RandomMutsVsCosineSimMLdepicted800.png")
Apoclass<-ggplot(allrands, aes(nmutsgbs, cosimdiag, col = APOBECHyper, fill = APOBECHyper))+
  geom_point(alpha=0.01)+
  stat_summary(fun.data = mean_se, geom="ribbon",alpha=0.5,)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  #scale_colour_viridis_b(option = "magma", breaks = c(0,3000,15000,30000))+ geom_hline(yintercept=0.95)+
  ylab("Cosine similarity \n to WGS profile")+
  xlab("Number of random \n mutations sampled")+
  theme(axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill = guide_legend(title="APOBEC \nsignature status"))
ggsave(plot=Apoclass,"figs/RandomMutsVsCosineSimAPOBEC_class.png",width=10,height=5)
Apoclass


expo<-NULL


i<-1
while(i<10){
  rando<-NULL
  print(paste("Processing", i))
  rando<-data.table(mclapply(seq(10,800,by=10), randsample, mc.cores=12))
  oss<-mclapply(rando$V1, mutinfo,mc.cores =8)
  rm(rando)
  allrands<-bind_rows(oss)
  expo[[i]]<-allrands
  rm(oss)
  rm(allrands)
  i<-i+1
}
then<-Sys.time()

timenew<-then - now


bleurgh<-data.table(bind_rows(expo))
bleurgh<-bleurgh[!is.nan(cosimdiag)]
saveRDS(bleurgh, "RandomEnzInfo800_scaledown.RDS")
bleurgh<-readRDS("RandomEnzInfo800_scaledown.RDS")
#bleurgh<-bleurgh[sample.int(nrow(bleurgh), 10000),]
bleurgh$APOBECHyper<-dfprops$APOBECHyper[match(bleurgh$Sample, dfprops$Var2)]
bleurgh$APOBECrate<-dfprops$Freq[match(bleurgh$Sample, dfprops$Var2)]
x<-110

accept<-0.95
takeoff<-lapply(seq(10,800,by =10), function(x){
  df<-subset(bleurgh, nmutsgbs==x)
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
  binomAh<-binom.test(Ah_success, Ah_tests, accept, alternative = "greater")
  resAh<-data.frame(APOBEChyper="APOBEC_hyper",
                    nmutsgbs=x,
                    p.val=binomAh$p.value,
                    lwr_95=binomAh$conf.int[1],
                    upr_95=binomAh$conf.int[2],
                    estimate=as.numeric(binomAh$estimate),
                    tests=as.numeric(binomAh$parameter),
                    prop=as.numeric(Ah_prop))
  
  binomAt<-binom.test(x=At_success, n=At_tests,accept, alternative = "greater")
  resAt<-data.frame(APOBEChyper="APOBEC_typical",
                    nmutsgbs=x,
                    p.val=binomAt$p.value,
                    lwr_95=binomAt$conf.int[1],
                    upr_95=binomAt$conf.int[2],
                    estimate=as.numeric(binomAt$estimate),
                    tests=as.numeric(binomAt$parameter),
                    prop=as.numeric(At_prop))
  binomAll<-binom.test(x=All_success, n=All_tests,accept, alternative = "greater")
  resAll<-data.frame(APOBEChyper="All_samples",
                     nmutsgbs=x,
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
saveRDS(takeoff,"takeoff.RDS")
takeoff$adj.pval<-p.adjust(takeoff$p.val)
first_sig_Ah<-takeoff$nmutsgbs[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_hyper")][1]
first_sig_Ah
first_sig_Ah_pval<-takeoff$adj.pval[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_hyper")][1]
first_sig_Ah_pval
takeoff$lwr_95[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_hyper")][1]
takeoff$upr_95[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_hyper")][1]
first_sig_At<-takeoff$nmutsgbs[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_typical")][1]
first_sig_At
first_sig_At_pval<-takeoff$adj.pval[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_typical")][1]
first_sig_At_pval
takeoff$lwr_95[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_typical")][1]
takeoff$upr_95[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "APOBEC_typical")][1]

first_sig_All<-takeoff$nmutsgbs[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "All_samples")][1]
first_sig_All
first_sig_All_pval<-takeoff$adj.pval[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "All_samples")][1]
first_sig_All_pval
takeoff$lwr_95[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "All_samples")][1]
takeoff$upr_95[which(takeoff$adj.pval<0.05 & takeoff$APOBEChyper == "All_samples")][1]

takeoff$APOBEChyper<-fct_relevel(takeoff$APOBEChyper, c("APOBEC_typical","APOBEC_hyper", "All_samples"))

takeoff_APO<-subset(takeoff, APOBEChyper != "All_samples")
prop_plot_APO<-ggplot(takeoff_APO, aes(nmutsgbs, prop, colour=APOBEChyper,fill=APOBEChyper)) + geom_bar(stat="identity", position="dodge")+#geom_point()+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  #geom_bar(stat="identity", position ="dodge")+
  geom_text(aes(x=first_sig_Ah,y=1.02, colour =APOBEChyper[1], label="*"), size =10, show.legend = FALSE)+
  geom_text(aes(x=first_sig_At,y=1.02, colour =APOBEChyper[2], label="*"), size =10,show.legend = FALSE)+
  xlab("Number of random \n mutations sampled")+
  ylab("Proportion of samples \nresulting in >90% \ncosine similarity") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"))
prop_plot_APO
ggsave(plot=prop_plot_APO,"figs/takeoff800APOBEC_combined.png",width=10,height=5)

takeoff_all<-subset(takeoff, APOBEChyper == "All_samples")
prop_plot_all<-ggplot(takeoff_all, aes(nmutsgbs, prop)) + geom_bar(stat="identity")+#geom_point()+
  #geom_bar(stat="identity", position ="dodge")+
  geom_text(aes(x=first_sig_All,y=1.02, label="*"), size =10,show.legend = FALSE)+#+
  xlab("Number of random \n mutations sampled")+
  ylab("Proportion of samples \nresulting in >90% \ncosine similarity") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(#plot.margin = unit(c(1.5,1,1,1), "lines"),
    axis.title = element_text(size = 18), axis.text =  element_text(size = 15),
    legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  guides(colour=guide_legend(title="APOBEC \nsignature status"),
         fill=guide_legend(title="APOBEC \nsignature status"))
prop_plot_all
ggsave(plot=prop_plot_all,"figs/takeoff800all_combined.png",width=10,height=5)
Randplotgrid<-Randplot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
grid1<-plot_grid(Randplotgrid, prop_plot_all, align = "v", axis = "lr",ncol=1)
ggsave(plot = grid1,"figs/Figure1.png", width =10, height =8)
grid1
Apoclassgrid<-Apoclass+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
grid2<-plot_grid(Apoclassgrid, prop_plot_APO,align="v",axis = "lr",ncol=1)
grid2
ggsave(plot = grid2,"figs/Figure2.png", width =10, height =8)
grid3<-plot_grid(MLplot,TMB_APO,align="v",axis = "lr",ncol=1)
grid3
ggsave(plot = grid3,"figs/Figure3.png", width =10, height =8)

