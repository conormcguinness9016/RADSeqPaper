# alldepth<-read.delim("ALLcov.txt", header = FALSE)
alldepth<-alldepth[alldepth$V3 > 50,]
alldepth10<-alldepth[alldepth$V1 == 10,]
alldep10lil<-A1A2depth10[which(alldepth10$V2>1120000 & alldepth10$V2<1220000),]
library(tidyverse)
ggplot(alldep10lil, aes(V2, V3)) + geom_line()
ggplot(alldepth10, aes(V2, V3)) + geom_line()
