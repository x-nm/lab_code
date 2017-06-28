# clean a GO download
# 2017.04.13 by xnm

setwd("D:\\xnm\\work\\17.3.16")

dat<-read.table("GO_ECM.txt",header = F)
dat<-apply(dat,1,toupper)
# dat<-as.data.frame(dat)
dup<-duplicated(dat)
dat<-subset(dat,!dup)
num<-length(dat)
dat<-c("GO_ECM",paste(">",num,sep = ""),dat)
write(dat,"GO_ECM.gmx",sep = "\t")
