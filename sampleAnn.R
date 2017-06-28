# 2017.2.20 by xnm

setwd("D:\\xnm\\work\\17.2.15\\KeyGene")
options(stringsAsFactors = F)

fileNm<-"WHP_NORM.txt"
# fileNm<-"WHP_NORM.txt"

dat<-read.table(fileNm, header = T)

# sample annotation file, for 
spName<-read.table("sampleAnn.txt",header = F)
index<-match(colnames(dat),spName[,1])
location<-!is.na(index)
for (i in 1:length(index)){
  if (location[i]){
    colnames(dat)[i]<-spName[index[i],2]
  }
}

write.table(dat,fileNm, quote = F, sep = "\t"， row.names = F)
