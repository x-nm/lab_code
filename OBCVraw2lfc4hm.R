# * in LFC FILE, ID NA were excluded.
# 2017.5.25 by xnm

setwd("D:\\xnm\\work\\DATA\\OB_CV\\raw")

########################
label<-"OB12"
########################

f<-paste0(label,"lfc_raw")
dat<-read.csv(f, header = T,sep = "\t", row.names = 1)

dat<-dat[,c(6,5,1)]
colnames(dat)<-c("ID","LFC","padj")

sum(dat$ID=="")
excluded<-which(dat$ID=="")
dat<-dat[-excluded,]

o<-paste0("../LFC_new/",label,"_LFC.txt")

write.table(dat,o,sep = "\t", quote = F)
