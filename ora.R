# 2017.3.9 by xnm
setwd("D:\\xnm\\work\\17.2.15\\KeyGene\\ORA/")
genesetFile<-"geneset_all.txt"
labelFile<-"label_all.txt"
# mergeLb<-"ORA" # Saving directory of ORA

# 3. over-representation analysis
# Input: genesetData(gsItem), sampleLabel(label[i])
# *gsItem is a vector, with geneID only
# *note: lfc file prepared is a file with full header (full 3 column names)
# *with the 1st col as gene ID. and cannot be row.names
# Output: "DEG in geneset","overlap#", "gsSize", "DEGsize","bgSize","pval"
ora<-function(gsData, spLabel){
  lfc<-read.table(paste(spLabel,"_LFC.txt",sep=""),header = T) # LFC
  DEG<-subset(lfc,padj<0.05)[,1] # DEG list
  DEGsize<-length(DEG)
  bgSize<-dim(lfc)[1]-DEGsize # NOTE, it's "number of other balls"
  gsSize<-length(gsData)
  DEGinSet<-DEG[na.omit(match(gsData,DEG))]
  overlap<-length(DEGinSet) # num of overlap
  DEGinSet<-paste(DEGinSet, collapse = ", ") # turn an array into a string
  pval<-phyper(q = overlap, 
               m = DEGsize , # number of red balls
               n = bgSize, # number of other balls
               k = gsSize, # number of balls drawn
               lower.tail = F) # hypergeometric test, SIGMA x>q
  # to calculate the statistical significance of having drawn a specific
  # k successes (out of n total draws) from the aforementioned population.
  res<-c(DEGinSet,overlap,gsSize,DEGsize,bgSize,pval)
  return(res)
}


# mergeLFCdata<-read.table(paste(mergeLb,"_LFC.txt",sep = ""),header = T)
# mergeNMdata<-read.table(paste(mergeLb,"_NORM.txt", sep = ""),header = T)

geneset<-readLines(genesetFile)
gsNum<-length(geneset)
label<-readLines(labelFile)
lbNum<-length(label)

# 2. for each sample, make an ora table
for (j in 1:lbNum){
  sample<-label[j]
  ORA<-c()
  
  # for each geneset, make a row of ora table
  for (i in 1:gsNum){
    gsItem<-readLines(paste("../geneset/",geneset[i],".gmx",sep = ""))
    dat<-c(gsItem[1],gsItem[2],ora(gsItem[-c(1,2)], sample)) 
    ORA<-rbind(ORA,dat)
  }
  colnames(ORA)<-c("GeneSet","Description","overlap","count of overlap","size of geneset", "size of DEG","size of background","pval" )
  ORA<-as.data.frame(ORA)
  p<-as.vector(ORA$pval)
  padj<-p.adjust(p, method = "fdr") # p is an array of p-vals
  ORA<-cbind(ORA,padj)
  
  ORA<-ORA[order(ORA$padj),]
  write.table(ORA, paste(sample,"_ORA.txt",sep = ""), quote = F, sep = "\t", row.names = F)
}


for (j in 1:lbNum){
  dat<-read.table(paste("ora/",label[j],"_ORA.txt",sep = ""),header = T)
}
