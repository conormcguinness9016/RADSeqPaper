setwd("./depth/")
library(stringr)
depth<-data.frame(sample=NA, depth=NA, stdev=NA)
sample <- NULL
depth <-NULL
stdev<-NULL
greater20<-NULL
# A1A12<-read.delim(list.files()[1])
# length(which(A1A12$X1 > 20))
list.files()
for(i in seq_along(list.files())){
  sample[i]<-str_extract(list.files()[i],"\\w*")
  print(sample[i])
  file<-read.delim(list.files()[i], header =FALSE)
  depth[i]<-median(file$V3)
  stdev[i]<-sd(file$V3)
  greater20[i]<-length(which(file$V3>20))
}
df<-data.frame(sample, depth, stdev, greater20)
df
write.table(df, "../depthtab.txt", sep = "\t", row.names = FALSE)
setwd("..")
A1A2depth<-read.delim("A1A2.bam.txt", header = FALSE)

A1A2depth10<-A1A2depth[A1A2depth$V1 == 10,]
A1A2depth10<-A1A2depth10[A1A2depth10$V3 > 20,]
A1A2dep10lil<-A1A2depth10[which(A1A2depth10$V2>1120000 & A1A2depth10$V2<1220000),]
plot(A1A2depth10$V2, log(A1A2depth10$V3), alpha = 0.1)
ggplot(A1A2depth10, aes(V2, V3)) + geom_line()
ggplot(A1A2dep10lil, aes(V2, V3)) + geom_point()
rm(A1A12depth)
q("yes")
