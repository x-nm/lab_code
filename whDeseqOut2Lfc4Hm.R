# WH DESEQ2 OUTPUT TO LFC
# 2017.5.25 by xnm

setwd("D:\\xnm\\work\\DATA\\WH\\deseq2_output")

############################
label<-"HN_Aorta"
############################

f<-paste0(label,"_LFC.txt")
dat<-read.csv(f,header = T, sep = "\t")

dat<-dat[,c(1,3,7)]
colnames(dat)<-c("ID","LFC","padj")

o<-paste0("LFC_new/",label,"_LFC.txt")
write.table(dat, o, sep = "\t", quote = F)



# # for PO
# setwd("D:\\xnm\\work\\DATA\\PO")
# dat<-read.csv("PO_LFC_old.txt", header = T, sep = "\t")
# colnames(dat)<-c("ID","LFC","padj")
# write.table(dat, "PO_LFC.txt", sep = "\t", quote = F)
