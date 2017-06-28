# usage: Rscript KeyGene.R geneset.txt label.txt
# *note: lfc file prepared is a file with full header (full 3 column names)
# *so as all the data matrix files; when read, with header = T, 
# *the 1st col is NAMED "ID" 
# 2017.2.16 by xnm



# REFINE: with "# for refining"
# merge data条目数太少，只有14621个；
# 看看每次match时na的数目，
# 准备一下全集的mergedata
# !!:na=gray heatmap
# !! match na也可能来自于基因名的不对应性！


# 1. USE WGCNA, NO CLUSTERING
# library("WGCNA")
# labeledHeatmap(Matrix = mat,
#                xLabels = colnames(mat),
#                yLabels = row.names(mat),
#                ySymbols = row.names(mat),
#                colorLabels = T,
#                colors = blueWhiteRed(100),
#                setStdMargins = F,
#                cex.lab.y = 0.9,
#                cex.lab.x = 0.9,
#                zlim = c(0,20),
#                main = paste("Heatmap of expression value of ",label, sep = ""))

# 2. pheatmap, with clustering, and pretty good picture 
# library(pheatmap)



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#                                DEFINATION                                   #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# setwd("D:\\Files\\Study\\交大\\课题\\17.2.15\\KeyGene")
setwd("D:\\xnm\\work\\17.2.15\\KeyGene")
genesetFile<-"geneset.txt"
labelFile<-"label.txt"
mergeLb<-"WHPA" # input data matrix is WHP_NORM.txt, WHP_LFC.txt

library(pheatmap)

# make a new dir named mergeLb
if (!dir.exists(mergeLb)){
  dir.create(mergeLb)
}

# Args<-commandArgs(trailingOnly = T)
# genesetFile<-Args[1]
# labelFile<-Args[2]

# Merge dataframe for all samples
# 怎么按sample聚类？或者，先不按样本聚类？

# TEST
# geneset<-readLines(genesetFile)
# gsNum<-length(geneset)
# i<-1
# gsData<-readLines(paste("geneset/",geneset[i],".gmx",sep = ""))
# gsData<-gsData[-c(1,2)]
# exprData<-read.table(paste(mergeLb,"_NORM.txt", sep = ""),header = T)
# lfcData<-read.table(paste(mergeLb,"_LFC.txt",sep = ""),header = T)
# label<-"test"



####################################### Functions #############################
# 1. Heatmap of gene expr (log(normValue+1))
# Input: genesetData(gsItem), exprDataframe, label for naming
# *gsItem is a vector, with geneID only
# Output: heatmap with cluster of genes and cluster of ALL SAMPLES for ONE GENESET
# *return match na for refining
exprHeatmap<-function(gsData, exprData, label){
  # get index
  index<-match(gsData, exprData$ID)
  na<-sum(is.na(index)) # for refining
  index<-na.omit(index)
  
  # remake the matrix
  mat<-exprData[index,]
  row.names(mat)<-mat$ID
  mat<-mat[,-1]
  # log(normValue+1)
  mat<-as.matrix(mat)
  mat<-log2(mat+1)
  
  # ?? set na as 0... will be na in merge all=T, cannot do the clustering
  # better to find na=gray
  mat[is.na(mat)]<-0

  # scale painting table
  wid<-440
  scale<- nrow(mat)/ncol(mat)
  if(scale>2.5){
    ht<-wid*scale*0.5
  }
  else if(scale<1){
    ht<-wid*scale*1.6
  }
  else{
    ht<-wid*scale
  }
  title <- paste(label,"_EHM", sep = "")
  bmp(paste(mergeLb,"/",title,".bmp",sep=""), width = wid, height = ht)
  
  # plotting
  # pheatmap(mat,fontsize = 5) # for plotting in studio
  pheatmap(mat, border_color = NA)
  
  dev.off()
  
  return(na) # for refining
}

# 2. Heatmap of lfc
# Input: genesetData(gsItem), lfcDataframe, label for naming
# *gsItem is a vector, with geneID only
# *lfc data: ID,lfc,padj
# Output: heatmap with cluster of genes and cluster of ALL SAMPLES for ONE GENESET
lfcHeatmap<-function(gsData, lfcData, label){
  # get index
  index<-na.omit(match(gsData, lfcData$ID))
  # remake the matrix
  mat<-lfcData[index,-c(3,5,7)] # exclude the padj col
  row.names(mat)<-mat$ID
  mat<-mat[,-1]
  mat<-as.matrix(mat)
  
  # ?? set na as 0... will be na in merge all=T, cannot do the clustering
  mat[is.na(mat)]<-0
  
  # scale painting table
  wid<-440 
  scale<- nrow(mat)/ncol(mat)
  ht<-wid*2  
  title <- paste(label,"_LHM", sep = "")
  bmp(paste(mergeLb,"/",title,".bmp",sep=""), width = wid, height = ht)
  # plotting
  pheatmap(mat, border_color = NA)
  dev.off()
}

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
###############################################################################



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#                                   Main                                      #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
mergeLFCdata<-read.table(paste(mergeLb,"_LFC.txt",sep = ""),header = T)
mergeNMdata<-read.table(paste(mergeLb,"_NORM.txt", sep = ""),header = T)

geneset<-readLines(genesetFile)
gsNum<-length(geneset)
label<-readLines(labelFile)
lbNum<-length(label)

# 1. for each geneset, plot 2 heatmaps
mna<-c() # for refining
for (i in 1:gsNum){
  gsLabel<-geneset[i]
  gsItem<-readLines(paste("geneset/",gsLabel,".gmx",sep = ""))
  gsItem<-gsItem[-c(1,2)] # exclude name and annotation of the geneSet
  hmLabel<-paste(mergeLb,"_",gsLabel,sep="")
  # plotting
  matchNa<-exprHeatmap(gsItem, mergeNMdata,hmLabel) # for refining
  lfcHeatmap(gsItem, mergeLFCdata,hmLabel)
  
  # for refining
  mna<-rbind(mna,c(gsLabel,matchNa,length(gsItem)))
}
# for refining
colnames(mna)<-c("geneset","na #","all #")
write.table(mna,paste(mergeLb,"/matchNa",sep=""),quote = F, sep = "\t") 




# 2. for each sample, make an ora table
for (j in 1:lbNum){
  sample<-label[j]
  ORA<-c()
  
  # for each geneset, make a row of ora table
  for (i in 1:gsNum){
    gsItem<-readLines(paste("geneset/",geneset[i],".gmx",sep = ""))
    dat<-c(gsItem[1],gsItem[2],ora(gsItem[-c(1,2)], sample)) 
    ORA<-rbind(ORA,dat)
  }
  colnames(ORA)<-c("GeneSet","Description","overlap","count of overlap","size of geneset", "size of DEG","size of background","pval" )
  ORA<-as.data.frame(ORA)
  p<-as.vector(ORA$pval)
  padj<-p.adjust(p, method = "fdr") # p is an array of p-vals
  ORA<-cbind(ORA,padj)
  
  ORA<-ORA[order(ORA$padj),]
  write.table(ORA, paste(mergeLb,"/",sample,"_ORA.txt",sep = ""), quote = F, sep = "\t", row.names = F)
}

