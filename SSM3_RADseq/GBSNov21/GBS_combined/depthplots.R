library(tidyverse)
chroms<-c(seq(1,19),"X","Y")
for(i in chroms){
  df<-read.delim(paste0("depthlong/depth_",i,".txt"),header=FALSE)
  names(df)<-c("sample","chrom","pos","depth")
  g<-ggplot(df, aes(pos, log(depth), col=sample))+geom_line(alpha=0.5)+facet_wrap(sample~.)+  
    ggtitle(paste0("Chromosome", i))+xlab("Position")+ylab("Depth (log scale)")+
    theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
          axis.title = element_text(size=25),
          legend.text = element_text(size=20),legend.title = element_text(size=22), strip.text = element_text(size=15))
  ggsave(plot=g, paste0("depthlong/",i,"_depthplot.pdf"),width=20,height=20)
}
  