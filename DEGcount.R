# count DEG
# 2017.5.25 by xnm



#################################
fdr<-0.1
fc<-0.5
#################################

setwd("D:\\xnm\\work\\17.5.25\\LFC")
spFiles<-readLines("samples.txt")

#################################
# DESEQ2
degCount<-function(f,fdr,fc)
{
  dat<-read.csv(paste0(f,"_LFC.txt"), header = T, sep = "\t")
  dat<-subset(dat,dat$padj<fdr)
  dat<-subset(dat,dat$LFC>log2(fc) | dat$LFC < log2((1/fc)))
  return(nrow(dat))
}

degCt<-c()
for (i in 1:length(spFiles)){
  deg<-degCount(spFiles[i],fdr,fc)
  degCt<-c(degCt,deg)
}
degCt<-cbind(spFiles,degCt)
colnames(degCt)<-c("samples","degCount")
write.table(degCt,paste0("../genesetSel/DEGCt_fdr",fdr,"fc",fc,".txt"),sep = "\t", quote = F, row.names = F)


#################################

#################################
# # DESEQ
# setwd("D:\\xnm\\work\\Rabbit_RNASeq\\DESeq_all")
# 
# degCount<-function(f){
#   dat<-read.csv(f, header = T, sep = "\t")
#   dat<-subset(dat,dat$padj<fdr)
#   dat<-subset(dat,dat$foldChange>fc | dat$foldChange < (1/fc))
#   print(f)
#   print(dim(dat)[1])
# }
# 
# # foldChange, padj
# files<-list.files()
# for (i in 1:length(files)){
#   degCount(files[i])
# }
#################################

