library(reticulate)
library(SigProfilerExtractorR)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerPlottingR)
b<-SigProfilerMatrixGeneratorR("BRCA", "mm10", "/Volumes/archive/cancergeneticslab/ConorM/GBSNov21/GBS_combined/totalVCFhets/", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
plotSBS(b + "96", "sigprofiler/","BRCA", "96",FALSE)
b<-SigProfilerMatrixGeneratorR("BRCA", "mm10", "/Volumes/archive/cancergeneticslab/ConorM/GBSNov21/GBS_combined/uniqueVCFhets/", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
plotSBS(b + "96", "sigprofiler/","BRCA", "96",FALSE)
d<-read.delim("unSNPS.txt", sep = " ", header = FALSE)
library(stringr)
library(devtools)
library(tidyverse)
library(forcats)
library(devtools)
library(SomaticSignatures)
d$V8<-as.numeric(str_extract(d$V7, "0\\.\\d{1,2}"))

colnames(d)<-c("CHROM", "POS", "REF", "ALT", "sample", "DP", "FORMAT", "AF")
dir.create("Rplotsgood")
png("Rplotsgood/varhist.png")
hist(d$AF)
dev.off()
d$clone<-"Monoclonal"
d$clone[str_detect(d$sample, "_D")]<-"Biclonal"
dcontrols<-d[str_detect(d$sample, "^D"),]
ggplot(dcontrols, aes(clone, AF))+geom_jitter()+geom_boxplot()

d$group_code<-str_extract(d$sample, "^\\w{1,1}")
d$group[d$group_code == "A"]<-"APOBEC3A + UGI"
d$group[d$group_code == "B"]<-"APOBEC3B + UGI"
d$group[d$group_code == "D"]<-"GFP+mCherry"
d$group[d$group_code == "O"]<-"0D5 (parental)"


subset(dmuts)
nmuts<-table(d$sample)

depth<-read.delim("depthtab.txt")
depth<-depth[!str_detect(depth$sample, "_D"),]
depth<-depth[depth$depth>20,]


depth<-depth[match(names(nmuts), depth$sample),]
depth$nmuts<-nmuts
depth<-na.omit(depth)
depth$mutrate<-(depth$nmuts*10^6)/depth$greater20
depth$group_code<-str_extract(depth$sample, "^\\w{1,1}")
depth$group[depth$group_code == "A"]<-"APOBEC3A + UGI"
depth$group[depth$group_code == "B"]<-"APOBEC3B + UGI"
depth$group[depth$group_code == "D"]<-"GFP+mCherry"
depth$group[depth$group_code == "O"]<-"0D5 (parental)"

depth$group<-fct_relevel(depth$group,"0D5 (parental)","GFP+mCherry","APOBEC3A + UGI","APOBEC3B + UGI")
dir.create("Rplotsgood")
png("Rplotsgood/mutationrate.png", width =1500, height=1500)
ggplot(depth, aes(group, mutrate)) + geom_boxplot()+
  geom_point(aes(col = sample,size= 40)) +
  xlab("Construct Transfected")+ylab("Mutation Rate \n (mutations/Mb sequenced)")+
  guides(color = guide_legend(override.aes = list(size = 10) ) )+
  theme(axis.title = element_text(size = 50), axis.text= element_text(size = 35),
  axis.text.x=element_text(angle = 45, hjust=1),
  legend.text=element_text(size =45), legend.title = element_text(size=50))
dev.off()
summary(aov(mutrate ~ group, data = depth))
d$start<-d$POS
d$end<-d$POS
mygenome <- FaFile("/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa")
dA<-subset(d, sample == "A2H1")
dG<-GRanges(d)
d_vr <- VRanges(
  seqnames = seqnames(dG),
  ranges = ranges(dG),
  ref = dG$REF,
  alt = dG$ALT,
  sampleNames = dG$sample,
  AF = dG$AF,
  group =dG$group)
mut_context<-mutationContext(d_vr, mygenome)
Mc<-data.frame(mut_context)
Mc$Context<-"Other"
Mc$Context[str_detect(Mc$context, "^T") & Mc$alteration=="CT"]<-"T(C>T)N"
chisq.test(Mc$Context, Mc$sampleNames)
aov(AF~Context*sampleNames, Mc)
aov(Context~sampleNames, Mc)
ggplot(Mc, aes(sampleNames, fill=Context, group=interaction(sampleNames,group)))+geom_bar(position = "fill")

ggplot(Mc, aes(sampleNames, AF, group=group, col=Context))+geom_boxplot()+ geom_jitter()
library(ggstatsplot)
png("Rplotsgood/fancybarplot.png", width =1500, height =1500)
ggbarstats(Mc,
  Context,
   sampleNames)
dev.off()
