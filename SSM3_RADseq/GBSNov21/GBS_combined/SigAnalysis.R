
library(BSgenome.Mmusculus.UCSC.mm10)

library("biomaRt")
library(SomaticSignatures)
library(tidyverse)
library(reshape2)
library(factoextra)
library(stringr)
library(MutationalPatterns)
mygenome <- FaFile("/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa")
setwd("GBS_combined/")
SSM3genoms<-read.delim("AllSNPsnoheader.txt", sep =" ", header = FALSE)

names(SSM3genoms)<-c("chr", "start", "ref", "alt", "sample", "FMT")
SSM3genoms$freqs<-str_extract(SSM3genoms$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
SSM3genoms<-SSM3genoms[-6]
SSM3genoms$RefFreq<-as.numeric(str_extract(SSM3genoms$freqs, pattern = "\\w{1,}(?=,)"))
SSM3genoms$AltFreq<-as.numeric(str_extract(SSM3genoms$freqs, "(?<=,)\\w{1,}"))
SSM3genoms$cov<-SSM3genoms$RefFreq + SSM3genoms$AltFreq
SSM3genoms$VAF<-SSM3genoms$AltFreq/SSM3genoms$cov

dir.create("Rplots")
hist(SSM3genoms$VAF[SSM3genoms$VAF!=1])
hist(SSM3genoms$VAF[SSM3genoms$VAF<0.9])
SSM3genoms$end<-SSM3genoms$start
SSM3genoms$group<-NA
SSM3genoms$group[str_detect(SSM3genoms$sample, "^A")]<-"A3A"
SSM3genoms$group[str_detect(SSM3genoms$sample, "^B")]<-"A3B"
SSM3genoms$group[str_detect(SSM3genoms$sample, "^D")]<-"NT"
SSM3genoms$clone<-"mono"
SSM3genoms$clone[str_detect(SSM3genoms$sample, "_D")]<-"mixed"
SSM3genoms_ss<-subset(SSM3genoms, cov >19 & sample != "B2F11")
SSM3genoms_ss$sample<-as.character(SSM3genoms_ss$sample)
SSM3genoms.gr<-makeGRangesFromDataFrame(SSM3genoms_ss, keep.extra.columns = TRUE)
png("Rplots/sharedvarshist.png", width =1500, height =1500)
hist(table(SSM3genoms.gr))
dev.off()
##what is the distribution of shared variants?
##most sites are unique
##filter out sites found in greater 10 samples
##find sites present in 10 or more samples


SSM3genoms_ssfilt<-SSM3genoms_ss %>% group_by(chr, start) %>% filter(n()<10)
hist(SSM3genoms_ssfilt$VAF[SSM3genoms_ssfilt$VAF<0.9])

SSM3genoms.gr<-makeGRangesFromDataFrame(SSM3genoms_ssfilt, keep.extra.columns = TRUE)
prop.table(table(SSM3genoms_ss$sample, SSM3genoms_ss$VAF>0.9),1)

##get the unique sites

SSM3_vr <- VRanges(
  seqnames = seqnames(SSM3genoms.gr),
  ranges = ranges(SSM3genoms.gr),
  ref = SSM3genoms.gr$ref,
  alt = SSM3genoms.gr$alt,
  sampleNames = SSM3genoms.gr$sample,
  group = SSM3genoms.gr$group)
SSM3_motifs<-mutationContext(SSM3_vr, mygenome)
SSM3_mm = motifMatrix(SSM3_motifs, group = "sampleNames", normalize = TRUE)
dir.create("Rplots")

png("Rplots/96profilebysample.png",width = 1000, height = 1000)
#plot_96_profile(SSM3_mm)
plotMutationSpectrum(SSM3_motifs, "sampleNames")
dev.off()
png("Rplots/cos_heatmap_all.png", width = 1000, height = 1000)
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3_mm), data.matrix(SSM3_mm)))
dev.off()

SSM3genoms_ssfilt<-SSM3genoms_ss %>% group_by(chr, start) %>% filter(n()<2)

SSM3genomsunique.gr<-makeGRangesFromDataFrame(SSM3genoms_ssfilt, keep.extra.columns = TRUE)
table(SSM3genoms$sample)

SSM3unique_vr <- VRanges(
  seqnames = seqnames(SSM3genomsunique.gr),
  ranges = ranges(SSM3genomsunique.gr),
  ref = SSM3genomsunique.gr$ref,
  alt = SSM3genomsunique.gr$alt,
  sampleNames = SSM3genomsunique.gr$sample,
  group = SSM3genomsunique.gr$group)
SSM3unique_motifs<-mutationContext(SSM3unique_vr, mygenome)
SSM3unique_mm = motifMatrix(SSM3unique_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysampleunique.png", width = 1000, height = 1000)
#plot_96_profile(SSM3unique_mm)
plotMutationSpectrum(SSM3unique_motifs, "sampleNames")
dev.off()
png("Rplots/cos_heatmap_unique.png", width = 1000, height = 1000)
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3unique_mm), data.matrix(SSM3unique_mm)))
dev.off()
##all profiles look the same-why
png("Rplots/VAFdist.png")
hist(SSM3genoms.gr$VAF)
dev.off()
##huge amount of subclonal muts in each sample-get rid of them as we are only interested in clonal mutations
SSM3genomhets.gr<-subset(SSM3genoms.gr, VAF >0.2 & VAF<0.7)

SSM3hets_vr <- VRanges(
  seqnames = seqnames(SSM3genomhets.gr),
  ranges = ranges(SSM3genomhets.gr),
  ref = SSM3genomhets.gr$ref,
  alt = SSM3genomhets.gr$alt,
  sampleNames = SSM3genomhets.gr$sample,
  group = SSM3genomhets.gr$group)
SSM3hets_motifs<-mutationContext(SSM3hets_vr, mygenome)
SSM3hets_mm = motifMatrix(SSM3hets_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysamplehets.png", width = 1000, height = 1000)
# plot_96_profile(SSM3hets_mm)

plotMutationSpectrum(SSM3hets_motifs, "sampleNames")

dev.off()
png("Rplots/cos_heatmap_hets.png", width = 1000, height = 1000)
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3hets_mm), data.matrix(SSM3hets_mm)))
dev.off()
dd<-subset(SSM3genomhets.gr, group == "NT")

png("Rplots/VAFbysamplehets.png")

ggplot(data.frame(dd), aes(clone, VAF)) +
geom_jitter()+geom_violin(alpha =0.5) +theme_classic()
dev.off()


A3Aonly<-subset(SSM3hets_motifs, group== "A3A")
sampleNames(A3Aonly)<-as.character(sampleNames(A3Aonly))
SSM3hetsA3A_mm = motifMatrix(A3Aonly, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysamplehetsA3Aonly.png", width = 1000, height = 1000)
# plot_96_profile(SSM3hetsA3A_mm)
plotMutationSpectrum(A3Aonly, "sampleNames")
dev.off()

SSM3genoms_ssfilt<-SSM3genoms_ss %>% group_by(chr, start) %>% filter(n()<6)
SSM3genomhetsunique.gr<-makeGRangesFromDataFrame(SSM3genoms_ssfilt, keep.extra.columns = TRUE)
SSM3genomhetsunique.gr<-subset(SSM3genomhetsunique.gr, VAF>0.2 & VAF<0.7)

SSM3genomhetsunique.gr<-unique(SSM3genomhets.gr)

SSM3hetsunique_vr <- VRanges(
  seqnames = seqnames(SSM3genomhetsunique.gr),
  ranges = ranges(SSM3genomhetsunique.gr),
  ref = SSM3genomhetsunique.gr$ref,
  alt = SSM3genomhetsunique.gr$alt,
  sampleNames = SSM3genomhetsunique.gr$sample,
  group = SSM3genomhetsunique.gr$group)
SSM3hetsunique_motifs<-mutationContext(SSM3hetsunique_vr, mygenome)
SSM3hetsunique_mm = motifMatrix(SSM3hetsunique_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysamplehetsunique.png", width = 1000, height = 1000)
#plot_96_profile(SSM3hetsunique_mm)
plotMutationSpectrum(SSM3hetsunique_motifs, "sampleNames")
dev.off()

# as.vector(table(match(SSM3genomhets.gr, SSM3genomhetsunique.gr))) < 2
# s3gu.gr<-SSM3genomhetsunique.gr[]
# SSM3hetsunique_vr <- VRanges(
#   seqnames = seqnames(s3gu.gr),
#   ranges = ranges(s3gu.gr),
#   ref = s3gu.gr$ref,
#   alt = s3gu.gr$alt,
#   sampleNames = s3gu.gr$sample,
#   group = s3gu.gr$group)
# SSM3hetsunique_motifs<-mutationContext(SSM3hetsunique_vr, mygenome)
# SSM3hetsunique_mm = motifMatrix(SSM3hetsunique_motifs, group = "sampleNames", normalize = TRUE)
# png("Rplots/96profilebysamplehetsunique.png", width = 1500, height = 1000)
# plot_96_profile(SSM3hetsunique_mm)
# dev.off()
png("Rplots/cos_heatmap_hetsunique.png", width = 1500, height = 1000)
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3hetsunique_mm), data.matrix(SSM3hetsunique_mm)))
dev.off()

SSM3_mm_group = motifMatrix(SSM3hetsunique_motifs, group = "group", normalize = TRUE)
png("Rplots/96profilebygrouphetsunique.png", width = 1500, height = 1000)
#plot_96_profile(SSM3_mm_group)
plotMutationSpectrum(SSM3hetsunique_motifs,group ="group")

dev.off()
SSM3_mm_group = motifMatrix(SSM3hets_motifs, group = "group", normalize = TRUE)
png("Rplots/96profilebygrouphets.png", width = 1500, height = 1000)
plot_96_profile(SSM3_mm_group)
dev.off()
dd<-subset(SSM3genomhetsunique.gr, group == "NT")

png("Rplots/VAFbysamplehetsunique.png")
ggplot(data.frame(dd), aes(clone, VAF)) +
geom_jitter()+geom_violin(alpha =0.5) +theme_classic()
dev.off()

dd<-subset(SSM3genoms.gr, group == "NT")

png("Rplots/VAFbysample.png")
ggplot(data.frame(dd), aes(clone, VAF)) +
geom_jitter()+geom_violin(alpha =0.5) +theme_classic()
dev.off()


SSM3genomsc.gr<-subset(SSM3genoms.gr, VAF <0.3)

SSM3sc_vr <- VRanges(
  seqnames = seqnames(SSM3genomsc.gr),
  ranges = ranges(SSM3genomsc.gr),
  ref = SSM3genomsc.gr$ref,
  alt = SSM3genomsc.gr$alt,
  sampleNames = SSM3genomsc.gr$sample,
  group = SSM3genomsc.gr$group)
SSM3sc_motifs<-mutationContext(SSM3sc_vr, mygenome)
SSM3sc_mm = motifMatrix(SSM3sc_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysamplesc.png", width = 1500, height = 1000)
#plot_96_profile(SSM3sc_mm)
plotMutationSpectrum(SSM3sc_motifs, group = "sampleNames")
dev.off()

SSM3genomscunique.gr<-unique(SSM3genomsc.gr)
SSM3sc_vr <- VRanges(
  seqnames = seqnames(SSM3genomscunique.gr),
  ranges = ranges(SSM3genomscunique.gr),
  ref = SSM3genomscunique.gr$ref,
  alt = SSM3genomscunique.gr$alt,
  sampleNames = SSM3genomscunique.gr$sample,
  group = SSM3genomscunique.gr$group)
SSM3sc_motifs<-mutationContext(SSM3sc_vr, mygenome)
SSM3sc_mm = motifMatrix(SSM3sc_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysampleunique_sc.png", width = 1500, height = 1000)
plotMutationSpectrum(SSM3sc_motifs, group = "sampleNames")
dev.off()



#
#
#
#
#
#
# n_sigs=4
# #' set the number of signatures to 4, and plot the signatures for both NMF and PCA
# sigs_nmf = identifySignatures(SSM3_mm, n_sigs, nmfDecomposition)
# plotSignatures(sigs_nmf) + ggtitle("NMF barchart")
#
# plot_96_profile(SSM3_mm)
#
# EditGroups_mm = motifMatrix(SSM3_motifs, group = "group", normalize = TRUE)
# EditGroups_mm
# plot_96_profile(EditGroups_mm)
#
# ##
#
# SSM3uniquegenoms<-read.delim("uniqueSNPsnoheader.txt", sep =" ", header = FALSE)
#
# names(SSM3uniquegenoms)<-c("chr", "start", "ref", "alt", "sample", "FMT")
#
#
#
# SSM3uniquegenoms$end<-SSM3uniquegenoms$start
# SSM3uniquegenoms$group<-NA
# SSM3uniquegenoms$group[str_detect(SSM3uniquegenoms$sample, "^A")]<-"APOBEC'd"
# SSM3uniquegenoms$group[str_detect(SSM3uniquegenoms$sample, "^E")]<-"E"
# SSM3uniquegenoms<-SSM3uniquegenoms[which(SSM3uniquegenoms$sample != c("E1D12") & SSM3uniquegenoms$sample != c("A1A12")),]
# SSM3uniquegenoms$sample<-as.character(SSM3uniquegenoms$sample)
# SSM3uniquegenoms$freqs<-str_extract(SSM3uniquegenoms$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
# SSM3uniquegenoms<-SSM3uniquegenoms[-6]
# SSM3uniquegenoms$RefFreq<-as.numeric(str_extract(SSM3uniquegenoms$freqs, pattern = "\\w{1,}(?=,)"))
# SSM3uniquegenoms$AltFreq<-as.numeric(str_extract(SSM3uniquegenoms$freqs, "(?<=,)\\w{1,}"))
# SSM3uniquegenoms$cov<-SSM3uniquegenoms$RefFreq + SSM3uniquegenoms$AltFreq
# mean(SSM3uniquegenoms$cov)
# SSM3uniquegenoms$VAF<-SSM3uniquegenoms$AltFreq/SSM3uniquegenoms$cov
# SSM3uniquegenoms.gr<-makeGRangesFromDataFrame(SSM3uniquegenoms, keep.extra.columns = TRUE)
#
#
#
#
# SSM3unique_vr <- VRanges(
#   seqnames = seqnames(SSM3uniquegenoms.gr),
#   ranges = ranges(SSM3uniquegenoms.gr),
#   ref = SSM3uniquegenoms.gr$ref,
#   alt = SSM3uniquegenoms.gr$alt,
#   sampleNames = SSM3uniquegenoms.gr$sample,
#   group = SSM3uniquegenoms.gr$group,
#   VAF=SSM3uniquegenoms.gr$VAF)
# SSM3unique_motifs<-mutationContext(SSM3unique_vr, mygenome)
#
# SSM3unique_mm = motifMatrix(SSM3unique_motifs, group = "sampleNames", normalize = TRUE)
# plot_96_profile(SSM3unique_mm)
# n_sigs=8
# #' set the number of signatures to 4, and plot the signatures for both NMF and PCA
# sigs_nmf = identifySignatures(SSM3unique_mm, n_sigs, nmfDecomposition)
# plotSignatures(sigs_nmf) + ggtitle("NMF barchart")
# plotSamples(sigs_nmf)
cosmic<-read.csv("~/Desktop/COSMIC/sigProfiler_SBS_signatures_2018_03_28.csv", stringsAsFactors = FALSE)

cosmicsigs<-data.matrix(cosmic[,3:37])
plot_96_profile(cosmicsigs[c(2,3,17)])
png("Rplots/sigs_vs_cosmic.png", width = 1500, height = 1000)
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3hets_mm), cosmicsigs))
dev.off()
EditGroups_mm = motifMatrix(SSM3unique_motifs, group = "group", normalize = TRUE)
EditGroups_mm
plot_96_profile(EditGroups_mm)
depthinf<-read.delim("depthtab.txt", header = TRUE, sep = "\t")
#colnames(depthinf)<-c("sample","greater20")
subset(SSM3genomhets.gr, )
depthinf<-depthinf[depthinf$sample %in% names(table(SSM3genomhets.gr$sample)),]
nmuts<-table(SSM3genomhets.gr$sample)
#nmuts<-nmuts[!str_detect(names(nmuts), "_D")]
nmuts

depthinf$nSNVs<-nmuts[match(names(nmuts), depthinf$sample)]
depthinf$mutrate<-depthinf$nSNVs/depthinf$greater20*1000000
depthinf$mutrate

depthinf
depthinf$group[str_detect(depthinf$sample, "^A")]<-"A3A"
depthinf$group[str_detect(depthinf$sample, "^B")]<-"A3B"
depthinf$group[str_detect(depthinf$sample, "^D")]<-"control"
depthinf$mutrate
depthinf
png("Rplots/mutratesallSNVs.png", width = 1500, height = 1000)
ggplot(depthinf, aes(group, mutrate))  +
  geom_boxplot(aes(alpha=0.5)) +
  geom_jitter(aes(colour = sample)) +
  ylab("Mutation rate (muts/Mb, all mutations)")+
  theme(text=element_text(size=20))
dev.off()

hets_mm_sample_nn = motifMatrix(SSM3hets_motifs, group = "sampleNames", normalize = FALSE)

A3Btargets<-colSums(hets_mm_sample_nn[c(13:16,29:32,45:48),])
depthinf$A3Btargetshets<-A3Btargets
depthinf$A3Btargethetsrate<-depthinf$A3Btargetshets/depthinf$greater20*1000000
depthinf
png("Rplots/A3btargethetrate.png", width = 1500, height = 1000)
ggplot(depthinf, aes(group, A3Btargethetsrate)) +
  geom_jitter(aes(colour = sample)) +
  geom_boxplot(alpha = 0.1)+
  ylab("Rate of APOBEC muts in TCN context (muts/Mb \n normalised to number of bp sequenced)")+
  theme(text=element_text(size=20))
dev.off()


###split into clonals
# hist(SSM3uniquegenoms$cov[SSM3uniquegenoms$sample == "A1A2"], binwidth =0.5)
# png("Rplots/VAFuniquevars.png", width = 1500, height = 1000, width = 1500)
# ggplot(SSM3uniquegenoms[which(SSM3uniquegenoms$VAF!=1),], aes(VAF))+geom_density(aes(colour = sample), bw =0.05, alpha =0.05)+geom_density(bw=0.05) +
#   +geom_rect()
#   theme_classic()
# dev.off()
# ?geom_rect()
# png("Rplots/VAFuniquevarsnoHoms.png", width =150)
# ggplot(SSM3genoms_ss[SSM3uniquegenoms$VAF!=1,], aes(VAF)) + geom_density(bw=0.025) +theme_classic()
# dev.off()
#
# png("Rplots/VAFuniquevarsnoHoms.png")
# ggplot(SSM3uniquegenoms, aes(VAF, colour = group)) + geom_density(bw=0.025) +theme_classic()
# dev.off()
#
# png("Rplots/VAFbygroup.png")
# ggplot(SSM3uniquegenoms[SSM3uniquegenoms$VAF!=1,], aes(VAF, colour =group)) + geom_density(bw=0.025) +theme_classic()
# dev.off()
#
# png("Rplots/VAFbysample.png")
# ggplot(SSM3genoms_ss[SSM3genoms_ss$VAF!=1,], aes(VAF, colour =sample)) + geom_density(bw=0.025) +theme_classic()
# dev.off()
#
# hets<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF<0.65 & SSM3uniquegenoms$VAF>0.35)),]
# depthinf$nSNVshets<-table(hets$sample)
# depthinf$mutratehets<-as.numeric(depthinf$nSNVshets)*10^6/depthinf$greater20
# png("Rplots/mutratehets.png")
# ggplot(depthinf, aes(group, mutratehets)) + geom_jitter(aes(colour = sample)) +
#   geom_boxplot(alpha = 0.5) + ylab("Heterozygous mutation rate (muts/Mb)")
# dev.off()
# hist(hets$VAF)
# hetsgenoms.gr<-makeGRangesFromDataFrame(hets, keep.extra.columns = TRUE)
# hets_vr <- VRanges(
#   seqnames = seqnames(hetsgenoms.gr),
#   ranges = ranges(hetsgenoms.gr),
#   ref = hetsgenoms.gr$ref,
#   alt = hetsgenoms.gr$alt,
#   sampleNames = hetsgenoms.gr$sample,
#   group = hetsgenoms.gr$group)
# hets_motifs<-mutationContext(hets_vr, mygenome)
# hets_mm_sample = motifMatrix(hets_motifs, group = "sampleNames", normalize = TRUE)
# png("Rplots/96profilebysamplehetsonly.png")
# plot_96_profile(hets_mm_sample)
# dev.off()
#
#
#
# hets_mm = motifMatrix(hets_motifs, group = "group", normalize = TRUE)
#
# png("Rplots/96profilebygrouphetsonly.png")
# plot_96_profile(hets_mm)
# dev.off()
#
#
# hets_mm_sample_nn = motifMatrix(hets_motifs, group = "sampleNames", normalize = FALSE)
# rownames(hets_mm_sample_nn)
# A3Btargets<-colSums(hets_mm_sample_nn[c(13:16,29:32,45:48),])
# depthinf$A3Btargetshets<-A3Btargets
# depthinf$A3Btargethetsrate<-depthinf$A3Btargetshets/depthinf$greater20
# depthinf
# png("Rplots/A3btargethetrate.png")
# ggplot(depthinf, aes(group, A3Btargethetsrate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb \n normalised to number of bp sequenced)")
# plot_96_profile(hets_mm_sample_nn)
# dev.off()
#
# subclonal<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF<0.35)),]
# depthinf$nSNVssubclonal<-table(subclonal$sample)
# depthinf$mutratesubclonal<-as.numeric(depthinf$nSNVssubclonal)*10^6/depthinf$greater20
# png("Rplots/mutratesubclonal.png")
# ggplot(depthinf, aes(group, mutratesubclonal)) + geom_boxplot(alpha = 0.5)+ geom_jitter(aes(colour = sample)) +
#     ylab("subclonal mutation rate (muts/Mb)")
# dev.off()
# hist(subclonal$VAF)
# subclonalgenoms.gr<-makeGRangesFromDataFrame(subclonal, keep.extra.columns = TRUE)
# subclonal_vr <- VRanges(
#   seqnames = seqnames(subclonalgenoms.gr),
#   ranges = ranges(subclonalgenoms.gr),
#   ref = subclonalgenoms.gr$ref,
#   alt = subclonalgenoms.gr$alt,
#   sampleNames = subclonalgenoms.gr$sample,
#   group = subclonalgenoms.gr$group)
# subclonal_motifs<-mutationContext(subclonal_vr, mygenome)
#
# subclonal_mm_sample_nn = motifMatrix(subclonal_motifs, group = "sampleNames", normalize = FALSE)
# A3BtargetsSubclone<-colSums(subclonal_mm_sample_nn[c(29:32,45:48),])
# depthinf$A3BtargetsSubclone<-A3BtargetsSubclone
# depthinf$A3BtargetSubclonerate<-depthinf$A3BtargetsSubclone*1000000/depthinf$greater20
# png("Rplots/A3btargetsubclonerate.png")
#
# ggplot(depthinf, aes(group, A3BtargetSubclonerate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb)")
# dev.off()
# subclonal_mm = motifMatrix(subclonal_motifs, group ="group", normalize = FALSE)
# png("Rplots/96profilesubclone.png")
# plot_96_profile(subclonal_mm)
# dev.off()
#
#
# homs<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF>0.65)),]
# depthinf$nSNVshoms<-table(homs$sample)
# depthinf$mutratehoms<-as.numeric(depthinf$nSNVshoms)*10^6/depthinf$greater20
# png("Rplots/mutratehoms.png")
# ggplot(depthinf, aes(group, mutratehoms)) + geom_jitter(aes(colour = sample)) +
#   geom_boxplot(alpha = 0.5) + ylab("homs mutation rate (muts/Mb)")
# dev.off()
# hist(homs$VAF)
# homsgenoms.gr<-makeGRangesFromDataFrame(homs, keep.extra.columns = TRUE)
# homs_vr <- VRanges(
#   seqnames = seqnames(homsgenoms.gr),
#   ranges = ranges(homsgenoms.gr),
#   ref = homsgenoms.gr$ref,
#   alt = homsgenoms.gr$alt,
#   sampleNames = homsgenoms.gr$sample,
#   group = homsgenoms.gr$group)
# homs_motifs<-mutationContext(homs_vr, mygenome)
#
# homs_mm_sample_nn = motifMatrix(homs_motifs, group = "sampleNames", normalize = FALSE)
# A3Btargetshoms<-colSums(homs_mm_sample_nn[c(29:32,45:48),])
# depthinf$A3Btargetshoms<-A3Btargetshoms
# depthinf$A3Btargethomsrate<-depthinf$A3Btargetshoms/depthinf$greater20
# ggplot(depthinf, aes(group, A3Btargethomsrate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb)")
# homs_mm = motifMatrix(homs_motifs, group ="group", normalize = FALSE)
# plot_96_profile(homs_mm)
# plot_96_profile(hets_mm)
# plot_96_profile(subclonal_mm)
#
#
# allsub<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF<0.3)),]
# depthinf$nSNVsallsub<-table(allsub$sample)
# depthinf$mutrateallsub<-as.numeric(depthinf$nSNVsallsub)*10^6/depthinf$greater20
# png("Rplots/mutrateallsub.png")
# ggplot(depthinf, aes(group, mutrateallsub)) + geom_boxplot(alpha = 0.5)+geom_jitter(aes(colour = sample)) +
#   ylab("allsub mutation rate (muts/Mb)")
# dev.off()
# hist(allsub$VAF)
# ggplot(allsub,aes(VAF)) + geom_density(aes(colour = sample))
# allsubgenoms.gr<-makeGRangesFromDataFrame(allsub, keep.extra.columns = TRUE)
# allsub_vr <- VRanges(
#   seqnames = seqnames(allsubgenoms.gr),
#   ranges = ranges(allsubgenoms.gr),
#   ref = allsubgenoms.gr$ref,
#   alt = allsubgenoms.gr$alt,
#   sampleNames = allsubgenoms.gr$sample,
#   group = allsubgenoms.gr$group)
# allsub_motifs<-mutationContext(allsub_vr, mygenome)
# allsub_mm_sample_nn = motifMatrix(allsub_motifs, group = "sampleNames", normalize = FALSE)
# A3Btargetsallsub<-colSums(allsub_mm_sample_nn[c(29:32,45:48),])
# depthinf$A3Btargetsallsub<-A3Btargetsallsub
# depthinf$A3Btargetallsubrate<-depthinf$A3Btargetsallsub/depthinf$greater20
# ggplot(depthinf, aes(group, A3Btargetallsubrate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb)")
# allsub_mm = motifMatrix(allsub_motifs, group ="group", normalize = FALSE)
# plot_96_profile(allsub_mm)
# plot_96_profile(hets_mm)
# plot_96_profile(subclonal_mm)
# depthinf
# ###split by group
# SSM3groupgenoms<-read.delim("SNPsbygroupnoheader.txt", sep =" ", header = FALSE)
#
# names(SSM3groupgenoms)<-c("chr", "start", "ref", "alt", "sample", "FMT")
#
#
#
# SSM3groupgenoms$end<-SSM3groupgenoms$start
# SSM3groupgenoms$group<-NA
# SSM3groupgenoms$group[str_detect(SSM3groupgenoms$sample, "^A")]<-"APOBEC'd"
# SSM3groupgenoms$group[str_detect(SSM3groupgenoms$sample, "^E")]<-"E"
# SSM3groupgenoms<-SSM3groupgenoms[which(SSM3groupgenoms$sample != c("E1D12") & SSM3groupgenoms$sample != c("A1A12")),]
# SSM3groupgenoms$sample<-as.character(SSM3groupgenoms$sample)
# SSM3groupgenoms$freqs<-str_extract(SSM3groupgenoms$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
# SSM3groupgenoms<-SSM3groupgenoms[-6]
# SSM3groupgenoms$RefFreq<-as.numeric(str_extract(SSM3groupgenoms$freqs, pattern = "\\w{1,}(?=,)"))
# SSM3groupgenoms$AltFreq<-as.numeric(str_extract(SSM3groupgenoms$freqs, "(?<=,)\\w{1,}"))
# SSM3groupgenoms$cov<-SSM3groupgenoms$RefFreq + SSM3groupgenoms$AltFreq
# (SSM3groupgenoms$cov)
# SSM3groupgenoms$VAF<-SSM3groupgenoms$AltFreq/SSM3groupgenoms$cov
# hist(SSM3groupgenoms$VAF)
# SSM3groupgenoms.gr<-makeGRangesFromDataFrame(SSM3groupgenoms, keep.extra.columns = TRUE)
#
#
#
#
# SSM3group_vr <- VRanges(
#   seqnames = seqnames(SSM3groupgenoms.gr),
#   ranges = ranges(SSM3groupgenoms.gr),
#   ref = SSM3groupgenoms.gr$ref,
#   alt = SSM3groupgenoms.gr$alt,
#   sampleNames = SSM3groupgenoms.gr$sample,
#   group = SSM3groupgenoms.gr$group)
# SSM3group_motifs<-mutationContext(SSM3group_vr, mygenome)
# SSM3group_motifs
# SSM3group_mm = motifMatrix(SSM3group_motifs, group = "sampleNames", normalize = FALSE)
# plot_96_profile(SSM3group_mm)
# A3Btargets<-colSums(SSM3group_mm[c(29:32,45:48),])
# data.frame(A3Btargets, names(A3Btargets))
# SSM3groupgenoms$VAF
# n_sigs=8
# #' set the number of signatures to 4, and plot the signatures for both NMF and PCA
# sigs_nmf = identifySignatures(SSM3group_mm, n_sigs, nmfDecomposition)
# plotSignatures(sigs_nmf) + ggtitle("NMF barchart")
# plotSamples(sigs_nmf)
# cosmic<-read.csv("~/Desktop/COSMIC/sigProfiler_SBS_signatures_2018_03_28.csv", stringsAsFactors = FALSE)
#
# cosmicsigs<-data.matrix(cosmic[,3:37])
# plot_96_profile(cosmicsigs[c(2,3,17)])
# plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3group_mm), cosmicsigs))
# EditGroups_mm = motifMatrix(SSM3group_motifs, group = "group", normalize = TRUE)
# EditGroups_mm
# plot_96_profile(EditGroups_mm)
# depthinf<-read.delim("depthtab.txt")
# depthinf<-depthinf[depthinf$sample %in% names(table(SSM3groupgenoms$sample)),]
# nmuts<-table(SSM3groupgenoms$sample)
#
# depthinf$nSNVs<-nmuts[match(data.frame(nmuts)$Var1, depthinf$sample)]
# depthinf$mutrate<-depthinf$nSNVs/depthinf$greater20*1000000
# depthinf$mutrate
#
# depthinf
# depthinf$group[str_detect(depthinf$sample, "^A")]<-"APOBEC'd"
# depthinf$group[str_detect(depthinf$sample, "^E")]<-"control"
# depthinf$mutrate
# group[str_detect(SSM3groupgenoms$sample, "^A")]<-"APOBEC'd"
# png("Rplots/mutratesallSNVs.png")
# ggplot(depthinf, aes(group, mutrate))  + geom_jitter() + geom_boxplot(aes(alpha=0.5)) + ylab("Mutation rate (muts/Mb, all mutations)")
# dev.off()
#
# SuM<-data.frame(SSM3unique_motifs)
#
# ggplot(SuM[which(SuM$VAF<0.9),], aes(VAF)) + geom_density(aes(colour = sampleNames), bw = 0.05, alpha = 0.5) + geom_density(bw = 0.05) + theme_classic()
# A3Btargets<-SuM[which(SuM$alteration == "CT" | SuM$alteration == "CG"),]
# hist(SuM$VAF[SuM$sampleNames == "A2D3"])
# hist(SuM$VAF[SuM$sampleNames == "E1D7"])
# A3Btargets<-A3Btargets[str_detect(string = A3Btargets$context, pattern = "^T"),]
# A3Btargets<-A3Btargets[which(A3Btargets$VAF !=1),]
# A3Btargets
# table(A3Btargets$sampleNames)/depthinf$greater20
# ggplot(A3Btargets, aes(VAF)) + geom_density(aes(colour = sampleNames), stat = "count")
# A3Bspl<-split(A3Btargets, A3Btargets$sampleNames)
# hist(A3Bspl$A2D3$VAF, bw= 0.05)
# normfactor<-depthinf$greater20
# hist
#

SSM3indels<-read.delim("uniqueindelsnoheader.txt", sep =" ", header = FALSE)
names(SSM3indels)<-c("chr", "start", "ref", "alt", "sample", "FMT")
SSM3indels$end<-SSM3indels$start
SSM3indels$group<-NA
SSM3indels$group[str_detect(SSM3indels$sample, "^A")]<-"A3A"
SSM3indels$group[str_detect(SSM3indels$sample, "^B")]<-"A3B"
SSM3indels$group[str_detect(SSM3indels$sample, "^D")]<-"control"

SSM3indels$sample<-as.character(SSM3indels$sample)
SSM3indels$freqs<-str_extract(SSM3indels$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
SSM3indels<-SSM3indels[-6]
SSM3indels$RefFreq<-as.numeric(str_extract(SSM3indels$freqs, pattern = "\\w{1,}(?=,)"))
SSM3indels$AltFreq<-as.numeric(str_extract(SSM3indels$freqs, "(?<=,)\\w{1,}"))
SSM3indels$cov<-SSM3indels$RefFreq + SSM3indels$AltFreq
mean(SSM3indels$cov)
SSM3indels$VAF<-SSM3indels$AltFreq/SSM3indels$cov
png("Rplots/indelVAF.png", width = 1500, height = 1000)
hist(SSM3indels$VAF)
dev.off()
SSM3indels_ss<-subset(SSM3indels, cov >19 & VAF>0.2 & VAF<0.7)
SSM3indels.gr<-makeGRangesFromDataFrame(SSM3indels, keep.extra.columns = TRUE)
nindels<-table(SSM3indels_ss$sample)
depthinf$nindels<-nindels[match(data.frame(nindels)$Var1, depthinf$sample)]
depthinf$indelmutrate<-depthinf$nindels*1000000/depthinf$greater20
depthinf$indelmutrate
depthinf
nindelshet<-table(SSM3indels[SSM3indels$VAF<0.65 & SSM3indels$VAF> 0.35,]$sample)
depthinf$nindelshet<-nindelshet[match(data.frame(nindelshet)$Var1, depthinf$sample)]
depthinf$nindelshetmutrate<-depthinf$nindelshet*(10^6)/depthinf$greater20
png("Rplots/indelmutrate.png", width = 1500, height = 1000)
ggplot(depthinf, aes(group, indelmutrate, colour = sample))+ geom_boxplot() + geom_point()
dev.off()
