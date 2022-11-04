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
setwd("Nik_Zainal/")
print("loaded packages")
##scrapes html table of all available enzymes from NEB site
url<-getURL("https://international.neb.com/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes")

tables<-readHTMLTable(url, header = T)
print("scraped table")
##get just the restriction enzyme chartREl
RElist<-tables[[2]]
##cut out all the other crap in the chart and fix uo html table format
RElist<-RElist[c(1,3, 5:10)]
##add the actual column names
colnames(RElist)<-c("Enzyme", "Sequence", "Activity in Buffer 1.1", "Activity in Buffer 2.1", "Activity in Buffer 3.1", "Activity in CutSmart Buffer", "Inactivation Temp", "Incubation Temp")
##delete the column names in the middle of the chart
RElist<-RElist[!is.na(RElist[5]),]
RElist<-RElist[which(RElist$Enzyme!="Enzyme"),]
print("built df in R")
write.table(RElist,"initialEnz_no_filter.csv", sep=",")
##get image data for the methylation column
url<-htmlParse(getURL("https://international.neb.com/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes"))
##the path was extracted using Selector Gadget-user friendly CSS/XPath finder
path = '//td[(((count(preceding-sibling::*) + 1) = 14) and parent::*)]//img'
methdata<-xpathSApply(url, path, xmlAttrs)["src",]
##methdata only has 281 instances, while RElist has 285 rows, from the website there are 4 enzymes without methylation data- I-CeuI, I-SceI, PI-PspI, PI-SceI, will remove
RElist<-RElist[which(RElist$Enzyme!="I-CeuI" & RElist$Enzyme!= "I-SceI" &  RElist$Enzyme!="PI-PspI" & RElist$Enzyme!="PI-SceI"),]
print("excluded methylation enzymes")
RElist$methdata<-methdata
##kick out enzymes that are sensitive to methylation
RElist<-RElist[str_detect(methdata,"not-sensitive"),]
##let's also cut out enzymes that are only nicking (enzymes containing Nb/Nt at start)
RElist<-RElist[!str_detect(RElist$Enzyme, "Nb\\."),]
RElist<-RElist[!str_detect(RElist$Enzyme, "Nt\\."),]
Enzymes<-RElist$Sequence
names(Enzymes)<-RElist$Enzyme

##get just the enzyme name, remove high fidelity enzymes


names(Enzymes)<-str_extract(names(Enzymes), "\\w{1,}\\p{Uppercase}")
Enzymes<-Enzymes[!duplicated(names(Enzymes))]
####tst with just 10 enzymes
#Enzymes<-Enzymes[1:10]
source("Nikfraglibgenfunction.R")
enzdat<-enzyme_info_gen(Enzymes, enzymedata = TRUE)

enzdat
saveRDS(enzdat, "enzdat.RDS")
enzdat<-readRDS("enzdat.RDS")
write.table(enzdat,"enzdat.csv", sep = ",")
read.csv("enzdat.csv")

frags<-frag_lib_generation(enzdat,Hsapiens, chroms=c(1:24))
#saveRDS(frags, "frags.RDS")
frags<-readRDS("frags.RDS")




REfraglibs<-lapply(frags$REfraglibs, function(x){
  x$chr <- replace(x[,chr], which(x[,chr] == "chr24"), "chrY")
  x$chr <- replace(x[,chr], which(x[,chr] == "chr23"), "chrX")
  return(x)
})
rm(frags)

REfraglibs<-lapply(REfraglibs, function(x){
  x[chr %in% names(Hsapiens)[1:24],]
})
REfraglibs<-REfraglibs[lapply(REfraglibs, function(x){nrow(x)}) >0]

length(REfraglibs)
#saveRDS(REfraglibs, "REfraglibs.RDS")

saveRDS(REfraglibs, "REfraglibs.RDS")




mutinfofragsnew<-function(x, mutrate=FALSE, libs = NULL ){
  RE.gr=makeGRangesFromDataFrame(data.frame(x), keep.extra.columns = TRUE)
  overlap_motifs<-subsetByOverlaps(BRCA_motifs, RE.gr)
  overlap_mm = motifMatrix(overlap_motifs, group = "sampleNames", normalize = TRUE)
  overlapmatrix = motifMatrix(overlap_motifs, group = "sampleNames", normalize = FALSE)
  nmutswgs<-colSums(BRCAmatrix)
  df<-data.table(Sample = names(nmutswgs))
  df$nmutswgs<-nmutswgs
  l<-length(nmutswgs)
  om<-matrix(rep(0, nrow(overlap_mm)*length(nmutswgs)),  nrow=nrow(overlap_mm), ncol=length(nmutswgs))
  om<-data.table(om)
  btab<-as.data.table(overlap_motifs)
  btab$Context<-"Other"
  btab$Context[btab$context =="T.A" & btab$alteration=="CT"]<-"APOBEC"
  btab$Context[btab$context =="T.T" & btab$alteration=="CT"]<-"APOBEC"
  btab$Context[btab$context =="T.A" & btab$alteration=="CG"]<-"APOBEC"
  btab$Context[btab$context =="T.T" & btab$alteration=="CG"]<-"APOBEC"
  dfprops<-data.frame(prop.table(table(btab$Context, btab$sampleNames),2))
  dfprops<-subset(dfprops, Var1 == "APOBEC")
  dfprops$APOBECHyper[dfprops$Freq>0.4]<-"APOBEC_hyper"
  dfprops$APOBECHyper[dfprops$Freq<0.4]<-"APOBEC_typical"
  dfprops$APOBECHyper<-fct_relevel(dfprops$APOBECHyper,c("APOBEC_typical", "APOBEC_hyper"))
  overlap_mm<-data.table(overlap_mm)
  om[,which(names(nmutswgs) %in% colnames(overlap_mm))]<-overlap_mm
  overlap_mm<-om
  om<-matrix(rep(0, nrow(overlapmatrix)*length(nmutswgs)),  nrow=nrow(overlapmatrix), ncol=length(nmutswgs))
  om<-data.table(om)
  overlapmatrix<-data.table(overlapmatrix)
  om[,which(names(nmutswgs) %in% colnames(overlapmatrix))]<-overlapmatrix
  overlapmatrix<-om
  cosmat<-cos_sim_matrix(overlap_mm,BRCA_mm)
  df$cosimdiag<-c(diag(cosmat))
  df$cosimdiag[is.nan(df$cosimdiag)]<-0
  df$nmutsgbs<-colSums(overlapmatrix)
  if( typeof(x) == "list" && ncol(x)> 4){
    df$Library<-rep(as.character(x[1,5]), nrow(df))
  } else if( !is.null(libs) ) {
    df$Library<-rep(libs, nrow(df))
  } else {
    print("Warning: no library names provided")
  }
  if(mutrate == TRUE){
    gbscov<-sum(width(GenomicRanges::reduce(RE.gr)))
    wgscov<-sum(seqlengths(Hsapiens)[1:24])
    df$mutrategbs<-(df$nmutsgbs*10^6)/gbscov
    df$mutratewgs<-(nmutswgs*10^6)/wgscov
    df$APOBECHypergbs<-dfprops$APOBECHyper
    df$ptnamegbs<-dfprops$Var2
    df
  } else{
    print("Mutation rate not computed")
    df
  }
}





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


AllEnzInfo<-mclapply(REfraglibs, function(x){
  mutinfofragsnew(x, mutrate=TRUE)
}, mc.cores=10)



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

F1Info<-mutinfofragsnew(f1ranges, mutrate=TRUE, libs = "F1")
AllEnzInfo<-c(AllEnzInfo,list(F1Info))
names(AllEnzInfo)<-unlist(lapply(AllEnzInfo, function(x){x[1,5]}))



#saveRDS(AllEnzInfo, "AllEnzInfo.RDS")





