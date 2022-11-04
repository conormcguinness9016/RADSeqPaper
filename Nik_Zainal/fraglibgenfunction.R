
enzyme_info_gen<-function(enzymerecs, enzymedata=FALSE){
  df<-data.table(Sequence = enzymerecs)
  df$CutLocation[!str_detect(df$Sequence, "\\(")]<-"Within"
  df$CutLocation[str_detect(df$Sequence, "\\(")]<-"Outside"

  #now map where exactly it cuts-first for the "Within" enzymes
  df$FiveCutSite<-NA
  RECutSeqsWithin<-df$Sequence[df$CutLocation == "Within"]
  ReCutSeqSitesWithin<-data.table(str_locate(RECutSeqsWithin, "/"))
  df$FiveCutSite[df$CutLocation == "Within"]<-ReCutSeqSitesWithin$start-1
  df$ThreeCutSite<-df$FiveCutSite
  ##for those outside, first number in brackets indicates how far out the enzyme cuts from the end of the recognition sequence
  #need code to identify this-something like str_length + first number in brackets before "/" for the "Outside" cutters
  if( "Outside" %in% df$CutLocation){
    print("Outside cutter in")
    OutsideCutters<-df[df$CutLocation == "Outside",]
    OCseqs<-str_extract(OutsideCutters$Sequence, "\\w*")
    OCseqlengths<-str_length(OCseqs)
    OC5cuts<-str_extract(OutsideCutters$Sequence, "(?<=\\().{1,2}")
    OC5cuts<-as.numeric(str_replace(OC5cuts, "/", ""))
    OCFiveCutSites<-OCseqlengths+OC5cuts
    OC3cuts<-str_extract(OutsideCutters$Sequence, "(?<=\\/).{1,2}")
    OC3cuts<-as.numeric(str_replace(OC3cuts, "\\)", ""))
    OC3cuts
    OCThreeCutSites<-OCseqlengths+OC3cuts
    df[which(df$CutLocation=="Outside"),]$FiveCutSite<-OCFiveCutSites
    df[which(df$CutLocation=="Outside"),]$ThreeCutSite<-OCThreeCutSites
  } else {
    print("No outside cutter")
  }
  seqlist<-df$Sequence
  df$Seq<-str_replace(seqlist, "/", "")
  allseqlengths<-str_extract(df$Sequence, "\\w*")
  allseqlengths<-str_length(allseqlengths)
  seqs<-df$Seq
  nobracks<-str_extract(seqs, "\\w{1,}(?=\\()")
  nobracks<-nobracks[!is.na(nobracks)]
  seqs[str_detect(seqs, "\\w{1,}(?=\\()")]<-nobracks
  df$Seq<-seqs
  rm(seqs)
  df$fivesplit<-NA
  df$threesplit<-NA
  splitstrings<-str_split_fixed(seqlist, "/", 2)
  df[df$CutLocation=="Within",]$fivesplit<-splitstrings[df$CutLocation=="Within",1]
  df[df$CutLocation=="Within",]$threesplit<-splitstrings[df$CutLocation=="Within",2]
  ##for the Outside cutters can use above-paste N*number of bases away onto the sequence
  fiveprimeOC<-NULL
  threeprimeOC<-NULL
  if( "Outside" %in% df$CutLocation){
    for (i in 1:nrow(df[df$CutLocation=="Outside",])){
      if (OC5cuts[i]>0){
        Nter<-paste0(rep("N", OC5cuts[i]), collapse = "")
        fiveprimeOC[i]<-paste0(OCseqs[i], Nter)
        threeprimeOC[i]<-"N"
      }
      else {
        fiveprimeOC[i]<-str_sub(OCseqs[i], end = OC5cuts[i]-1)
        threeprimeOC[i]<-str_sub(OCseqs[i], start = OC5cuts[i])
      }
    }} else {
      print("I already told you, there aren't any outside cutters")
    }
  df$fivesplit[df$CutLocation=="Outside"]<-fiveprimeOC
  df$threesplit[df$CutLocation=="Outside"]<-threeprimeOC
  df$fivesplit[df$fivesplit==""]<-"N"
  g<-reverseComplement(DNAString(df$Seq[1]))
  df$revcompSeq<-as.character(reverseComplement(DNAStringSet(df$Seq)))
  print("All restriction enzyme data successfully loaded.")
  if(enzymedata == TRUE){
    df$Library<-names(enzymerecs)
  } else {
  df$Library<-df$Seq
}
  return(df)
}



frag_lib_generation<-function(df, genome, chroms=seq_along(genome), fragmin = 70, fragmax = 378){
  REfrags<-list()
  revcomp<-list()
  REfraglibs<-list()
  for (x in chroms) {
    #for (x in 21:22) {
    print(paste("Starting for chromosome", x))
    ##i defines different enzymes
    for(i in 1:nrow(df)){
      print(paste("Processing", df$Library[i], "on chromosome", x))
      ##match the forward pattern
      m<-Biostrings::matchPattern(df$Seq[i], genome[[x]], fixed = FALSE)
      ends<-end(gaps(m))
      ##cut where the enzyme cuts, not just the recognition site
      temp_df<-data.frame(end=ends + df$FiveCutSite[i] + 1,chr=paste0("chr", x))
      temp_df$start<-c(temp_df$end[1:nrow(temp_df)] - 1)
      print(paste("Processed forward strand!"))
      ##for enzymes with non palindromic recognition sequence, do the same with the reverse complement
      if(df$Seq[i] != df$revcompSeq[i]){
        m<-Biostrings::matchPattern(df$revcompSeq[i], genome[[x]], fixed = FALSE)
        ends<-end(gaps(m))
        revtemp_df<-data.frame(end=ends - df$ThreeCutSite[i] - 1, chr=paste0("chr", x))
        revtemp_df$start<-c(revtemp_df$end[1:nrow(revtemp_df)] - 1)
        revtemp_df<-revtemp_df[c("chr","start","end")]
        ##replace the ends with the ends of the chromosomes-problematic with human chromosomes as last ~1000bp is NN?
        revtemp_df$end<-replace(revtemp_df$end, which(revtemp_df$end>length(genome[[x]])), length(genome[[x]]))
        revtemp_df$start<- replace(revtemp_df$start, which(revtemp_df$start<0), 1)
        ##combine the forward and reverse recognition sites
        forandrev<-rbind(temp_df, revtemp_df)
        ##rearrange the data frame by the endpoints
        forandrev<-arrange(forandrev, end)
        forandrev$start<-c(1, forandrev$end[1:nrow(forandrev)-1] + 1)
        forandrev<-forandrev[c("chr","start","end")]
        forandrev$end<-replace(forandrev$end, which(forandrev$end>length(genome[[x]])), length(genome[[x]]))
        forandrev$start<- replace(forandrev$start, which(forandrev$start<0), 1)
        #find the sizes of the fragments
        forandrev$width<-forandrev$end-forandrev$start
        revcomp[[df$Library[i]]]<-forandrev
      }
      ##do the same for the palindromic recognition sites
      temp_df$start<-c(1, temp_df$end[1:nrow(temp_df)-1] + 1)
      temp_df<-temp_df[c("chr","start","end")]
      temp_df$end<-replace(temp_df$end, which(temp_df$end>length(genome[[x]])), length(genome[[x]]))
      temp_df$start<- replace(temp_df$start, which(temp_df$start<0), 1)
      temp_df$width<-temp_df$end-temp_df$start
      if(df$Seq[i] != df$revcompSeq[i]){
        ##add the data for this chromosome to the enzyme in the fragments list
        REfrags[[df$Library[i]]]<-data.table(rbind(forandrev, REfrags[[df$Library[i]]]))
      } else {
        REfrags[[df$Library[i]]]<-data.table(rbind(temp_df, REfrags[[df$Library[i]]]))
      }
      ##simulate "libraries" by restricting the size of the fragments
      REfraglibs[[df$Library[i]]]<-REfrags[[df$Library[i]]][which(REfrags[[df$Library[i]]][,4]>fragmin&REfrags[[df$Library[i]]][,4]<fragmax),]
      REfraglibs[[df$Library[i]]][,"Library"]<-rep(df$Library[i], nrow(REfraglibs[[df$Library[i]]]))
    }
    print(paste("Processed all Librarys for chromosome", x))

  }
  data.table(REfrags, REfraglibs)
}
