library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg19)

library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse)
length(Hsapiens)
library(SimRAD)
##not run-simulated DNA sequence-for developing tool
source("/Volumes/scratch/cancergeneticslab/ConorM/InSilicoRADseq/Rfunctions/fraglibgenfunction.R")
mousegenom<-readDNAStringSet("../../ref_genomes/mouse2/mouse.fa", format="fasta",
                nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

PstI<-"CTGCA/G"
MspI<-"C/CGG"
frags<-frag_lib_generation(enzyme_info_gen(c(PstI, MspI)), mousegenom)
mousegenom
now<-Sys.time()
chrs<-NULL
REfrags<-frags[[1]]

REfrags<-mapply(cbind, REfrags, "SampleID"=names(REfrags), SIMPLIFY=F)

AA<-rbind(REfrags[[1]], REfrags[[2]])
AA
AAo<-AA[order(AA[,1], AA[,3]),]
AAo
##update start coords
AAo$start<-c(1, AAo$end[1:nrow(AAo) - 1]+1)
AAo$width<-AAo$end-AAo$start
AAo
AAo<-AAo[which(AAo$width < 100 & AAo$width>39),]
nrow(AAo[which(AAo$SampleID == "CTGCAG"),])
library(IRanges)
library(Biostrings)

library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg19)

library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse)
length(Hsapiens)
library(SimRAD)

RElist$Enzyme<-as.character(RElist$Enzyme)
frag_lib_generation(RElist, mousegenom,fragmin = 37, fragmax = 100)











#simseq<-Hsapiens[[1]]
##Enzyme column needs to be changed from factor to character
RElist$Enzyme<-as.character(RElist$Enzyme)
##define empty lists for for loop

##fro x in 24: take chr1-22 and chrx and y
seq_along(mousegenom)
seq_along(mousegenom)


##replace chromosome 23+24 names with chrX and chrY
REfraglibs<-lapply(REfraglibs, function(x){
  x[,1] <- as.character(x[,1])
  x[,1] <- replace(x[,1], which(x[,1] == "chr24"), "chrY")
  x[,1] <- replace(x[,1], which(x[,1] == "chr23"), "chrX")
  return(x)
})
REfrags<-lapply(REfrags, function(x){
  x[,1] <- as.character(x[,1])
  x[,1] <- replace(x[,1], which(x[,1] == "chr24"), "chrY")
  x[,1] <- replace(x[,1], which(x[,1] == "chr23"), "chrX")
  return(x)
})
saveRDS(REfrags, "/Volumes/scratch/cancergeneticslab/ConorM/InSilicoRADseq/5merfrags.rds")
saveRDS(REfraglibs, "/Volumes/scratch/cancergeneticslab/ConorM/InSilicoRADseq/5merfraglibs.rds")
rm(REfrags)
rm(REfraglibs)
