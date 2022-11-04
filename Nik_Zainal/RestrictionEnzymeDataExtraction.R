library(XML)
library(RCurl)
library(stringr)
print("loaded packages")
##scrapes html table of all available enzymes from NEB site
url<-getURL("https://international.neb.com/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes")

tables<-readHTMLTable(url, header = T)
print("scraped table")
##get just the restriction enzyme chart
RElist<-tables[[2]]
##cut out all the other crap in the chart and fix uo html table format
RElist<-RElist[c(1,3, 5:10)]
##add the actual column names
colnames(RElist)<-c("Enzyme", "Sequence", "Activity in Buffer 1.1", "Activity in Buffer 2.1", "Activity in Buffer 3.1", "Activity in CutSmart Buffer", "Inactivation Temp", "Incubation Temp")
##delete the column names in the middle of the chart
RElist<-RElist[!is.na(RElist[5]),]
RElist<-RElist[which(RElist$Enzyme!="Enzyme"),]
print("built df in R")
##get image data for the methylation column
url<-htmlParse(getURL("https://international.neb.com/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes"))
##the path was extracted using Selector Gadget-user friendly CSS/XPath finder
path = '//td[(((count(preceding-sibling::*) + 1) = 14) and parent::*)]//img'
methdata<-xpathSApply(url, path, xmlAttrs)["src",]
##methdata only has 281 instances, while RElist has 285 rows, from the website there are 4 enzymes without methylation data- I-CeuI, I-SceI, PI-PspI, PI-SceI, will remove
RElist<-RElist[which(RElist$Enzyme!="I-CeuI" & RElist$Enzyme!= "I-SceI" &  RElist$Enzyme!="PI-PspI" & RElist$Enzyme!="PI-SceI"),]
print("excluded methylation enzymes")
##NEED TO FIX THIS UP TO BE READABLE
RElist$methdata<-methdata
##kick out enzymes that are sensitive to CpG methylation
RElist<-RElist[str_detect(methdata,"not-sensitive"),]
##let's also cut out enzymes that are only nicking (enzymes containing Nb/Nt at start)
RElist<-RElist[!str_detect(RElist$Enzyme, "Nb\\."),]
RElist<-RElist[!str_detect(RElist$Enzyme, "Nt\\."),]
##locate and log the cut sites of each enzyme
##need a way to split the enzymes between those that cut in the sequence and those that don't
##can do this as the Sequence column defines this through the presence of brackets
RElist$CutLocation<-NA
RElist[!str_detect(RElist$Sequence, "\\("),]$CutLocation<-"Within"
RElist[str_detect(RElist$Sequence, "\\("),]$CutLocation<-"Outside"
##now map where exactly it cuts-first for the "Within" enzymes
RElist$FiveCutSite<-NA
RECutSeqsWithin<-RElist[RElist$CutLocation == "Within",]$Sequence
ReCutSeqSitesWithin<-data.frame(str_locate(RECutSeqsWithin, "/"))
RECutSeqsWithin
RElist[RElist$CutLocation == "Within",]$FiveCutSite<-ReCutSeqSitesWithin$start-1
RElist$ThreeCutSite<-RElist$FiveCutSite
##for those outside, first number in brackets indicates how far out the enzyme cuts from the end of the recognition sequence
##need code to identify this-something like str_length + first number in brackets before "/" for the "Outside" cutters
OutsideCutters<-RElist[RElist$CutLocation == "Outside",]
OCseqs<-str_extract(OutsideCutters$Sequence, "\\w*")
OCseqlengths<-str_length(OCseqs)
OC5cuts<-str_extract(OutsideCutters$Sequence, "(?<=\\().{1,2}")
OC5cuts<-as.numeric(str_replace(OC5cuts, "/", ""))
OCFiveCutSites<-OCseqlengths+OC5cuts
OC3cuts<-str_extract(OutsideCutters$Sequence, "(?<=\\/).{1,2}")
OC3cuts<-as.numeric(str_replace(OC3cuts, "\\)", ""))
OC3cuts
OCThreeCutSites<-OCseqlengths+OC3cuts
RElist[which(RElist$CutLocation=="Outside"),]$FiveCutSite<-OCFiveCutSites
RElist[which(RElist$CutLocation=="Outside"),]$ThreeCutSite<-OCThreeCutSites
#hist(RElist$FiveCutSite[which(RElist$CutLocation == "Outside")])
#hist(RElist$FiveCutSite[which(RElist$CutLocation == "Within")])
library(tidyverse)

#ggplot(data=RElist, aes(x=FiveCutSite, fill=CutLocation)) + geom_density(alpha = 0.5) +xlab("Bases away from 5` end of recognition site")

##fix up the recognition sequences so they don't contain / or ()
library(stringr)
seqlist<-RElist$Sequence
RElist$Seq<-str_replace(seqlist, "/", "")
allseqlengths<-str_extract(RElist$Sequence, "\\w*")
allseqlengths<-str_length(allseqlengths)
seqs<-RElist$Seq
nobracks<-str_extract(seqs, "\\w{1,}(?=\\()")
nobracks<-nobracks[!is.na(nobracks)]
seqs[str_detect(seqs, "\\w{1,}(?=\\()")]<-nobracks
RElist$Seq<-seqs
rm(seqs)
RElist$Enzyme[str_detect(RElist$Seq, "TCA")]
RElist$Enzyme[str_detect(RElist$Seq, "AGT")]
RElist$Enzyme[str_detect(RElist$Seq, "AAT")]
##from the enzyme list, BclI would be good for sequencing an APOBECed genome (or BspCNI, BspHI)? Predict that more fragments would be sequenced if genome hadn't been APOBECed
##what about other signatures
##Signature 17
RElist$Enzyme[str_detect(RElist$Seq, "CTT")]


##finally we need to define the 5` and 3` breakpoints for each enzyme. For the "Within" enzymes this should be easy using str_splot
RElist$fivesplit<-NA
RElist$threesplit<-NA
splitstrings<-str_split_fixed(seqlist, "/", 2)
RElist[RElist$CutLocation=="Within",]$fivesplit<-splitstrings[RElist$CutLocation=="Within",1]
RElist[RElist$CutLocation=="Within",]$threesplit<-splitstrings[RElist$CutLocation=="Within",2]
##for the Outside cutters can use above-paste N*number of bases away onto the sequence
fiveprimeOC<-NULL
threeprimeOC<-NULL
for (i in 1:nrow(RElist[RElist$CutLocation=="Outside",])){
  if (OC5cuts[i]>0){
    Nter<-paste0(rep("N", OC5cuts[i]), collapse = "")
    fiveprimeOC[i]<-paste0(OCseqs[i], Nter)
    threeprimeOC[i]<-"N"
  }
  else {
    fiveprimeOC[i]<-str_sub(OCseqs[i], end = OC5cuts[i]-1)
    threeprimeOC[i]<-str_sub(OCseqs[i], start = OC5cuts[i])
  }
}
RElist[RElist$CutLocation=="Outside",]$fivesplit<-fiveprimeOC
RElist[RElist$CutLocation=="Outside",]$threesplit<-threeprimeOC
RElist$fivesplit[RElist$fivesplit==""]<-"N"

library(Biostrings)
g<-reverseComplement(DNAString(RElist$Seq[1]))
RElist$revcompSeq<-as.character(reverseComplement(DNAStringSet(RElist$Seq)))
print("All restriction enzyme data successfully loaded.")
