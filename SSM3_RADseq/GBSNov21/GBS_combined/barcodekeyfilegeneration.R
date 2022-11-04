####fast_x_barcode_splitter only accepts barcodes of the same length as input, so this script splits up the barcodes file into those of different lengths


library(stringr)
barcodes<-read.delim("key/translationtable.txt", header = TRUE)
bc<-data.frame(barcodes$uidtag,as.character(barcodes$barcode))
bc$barcodes.uidtag<-str_replace(bc$barcodes.uidtag, " ", "_")
seq<-c(0,rep(c(1,2,3), 15))
str
bc$fn<-paste0(bc$barcodes.uidtag, "_", seq)
cat(unlist(lapply(seq(nrow(bc)), function(x){
  paste0(">",bc[x,3], "\n","^",bc[x,2], "\n")
})), file = "bc.fasta")
write(unique(bc$barcodes.uidtag), file = "samples.txt")
