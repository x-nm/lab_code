# usage: Rscript KeyGene.R geneset.txt label.txt
# *note: lfc file prepared is a file with full header (full 3 column names)(ID, LFC, PADJ)
# *so as all the data matrix files; when read, with header = T, 
# *the 1st col is NAMED "ID" 
# *exprData is a matrix with ID as a col,
# *so in data preparation, the data.frame should be either with a plus col, or with a full header
# *Merged LFC data is with full col, and with padj col, too
# 2017.2.16 by xnm



# REFINE: with "# for refining"
# merge data条目数太少，只有14621个；
# 看看每次match时na的数目，
# 准备一下全集的mergedata
# !!:na=gray heatmap
# !! match na也可能来自于基因名的不对应性！


# for pheatmap dist NA problem:
# tryCatch():
# if catch error:
# draw 2 plots: one with NA<-0, the other with NA row (if NA included) removed
# <-better way: the other with dist NA row removed

# ！！如果2中只删掉Dup的那一半号码呢？
# 不过主要是对DLHM而言，毕竟DLHM主要看NA。那就先将就着看看DLHM1吧...因为NDEG set NA以后实在是太多空了..
# 不作hclust倒是可以看


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
genesetFile<-"geneset_all.txt"
labelFile<-"label.txt"
mergeLb<-"WHP" # input data matrix is WHP_NORM.txt, WHP_LFC.txt

library(pheatmap)
library(gplots)

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
# i<-2
# gsData<-readLines(paste("geneset/",geneset[i],".gmx",sep = ""))
# gsData<-gsData[-c(1,2)]
# exprData<-read.table(paste(mergeLb,"_NORM.txt", sep = ""),header = T)
# lfcData<-read.table(paste(mergeLb,"_LFC.txt",sep = ""),header = T)
# label<-"test"
# title<-"test"


####################################### Functions #############################
# 0.1 find the row index: rows with dist() val == NA
# input: the matrix
# output: matrix with those rows removed
matDelRows<-function(matr){
  distMt<-as.matrix(dist(matr))
  index<-which(is.na(distMt), arr.ind = T)
  index<-index[,1]
  index<-index[!duplicated(index)]
  matr<-matr[-index,]
  return(matr)
}

# 0.2 prepare a canvas and draw a heatmap
# Input: mat, title, width, height
hm<-function(matrix, title, width, height){
  depth<-round(max(abs(min(matrix, na.rm = T)),abs(max(matrix, na.rm = T))))
  c<-c(-(depth+2):(depth+1))
  len<-length(c)
  col<-bluered(len)
  
  bmp(paste(title,".bmp",sep=""), width = width, height = height)
  # plotting
  pheatmap(matrix, border_color = NA, breaks = c, color = col)
  dev.off()
  
  bmp(paste(title,"_num.bmp",sep=""), width = width, height = height)
  # plotting
  pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T) # TO SHOW THE NA
  dev.off()
}

ehm<-function(matrix, title, width, height){
  bmp(paste(title,".bmp",sep=""), width = width, height = height)
  # plotting
  pheatmap(matrix, border_color = NA)
  dev.off() 
}

# 1. Heatmap of gene expr (log(normValue+1))
# Input: genesetData(gsItem), exprDataframe, label for naming (saving dir included)
# *gsItem is a vector, with geneID only
# *exprData is a matrix with ID as a col
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

  # mat[is.na(mat)]<-0
  # delete all NA rows
  mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
  
  # draw the heatmap
  if(all(!is.na(dist(mat)))){
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
    ehm(mat, paste(label,"_EHM", sep = ""), wid, ht) # heatmap
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
    # scale painting table
    wid<-440
    scale<- nrow(mat1)/ncol(mat1)
    if(scale>2.5){
      ht<-wid*scale*0.5
    }
    else if(scale<1){
      ht<-wid*scale*1.6
    }
    else{
      ht<-wid*scale
    }
    ehm(mat1, paste(label,"_EHM1", sep = ""), wid, ht) # heatmap
    
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
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
    ehm(mat, paste(label,"_EHM2", sep = ""), wid, ht) # heatmap
  }

  return(na) # for refining
}

# 2. Heatmap of lfc
# Input: genesetData(gsItem), lfcDataframe, label for naming (saving dir included)
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

  # delete all NA rows
  mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
  
  # draw the heatmap
  if(all(!is.na(dist(mat)))){
    # scale painting table
    wid<-440 
    ht<-wid*2  
    hm(mat, paste(label,"_LHM", sep = ""), wid, ht) # heatmap
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
    # scale painting table
    wid<-440 
    ht<-wid*2  
    hm(mat1, paste(label,"_LHM1", sep = ""), wid, ht) # heatmap
    
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
    # scale painting table
    wid<-440 
    ht<-wid*2  
    hm(mat, paste(label,"_LHM2", sep = ""), wid, ht) # heatmap
  }
  
}

# 2.2 Heatmap of dlfc: with non-DEG colored as gray
dlfcHeatmap<-function(gsData, lfcData, label){
  # get index
  index<-na.omit(match(gsData, lfcData$ID))
  # remake the matrix
  mat<-lfcData[index,-c(3,5,7)]
  padj<-lfcData[index,-c(2,4,6)]
  mat[padj>0.05]<-NA
  # padj[is.na(padj)]<-0
  row.names(mat)<-mat$ID
  mat<-mat[,-1]
  mat<-as.matrix(mat)

  # delete all NA rows
  mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
  
  # draw the heatmap
  if(all(!is.na(dist(mat)))){
    # scale painting table
    wid<-440 
    ht<-wid*2  
    hm(mat, paste(label,"_DLHM", sep = ""), wid, ht) # heatmap
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
    # scale painting table
    wid<-440 
    ht<-wid*2  
    hm(mat1, paste(label,"_DLHM1", sep = ""), wid, ht) # heatmap
    
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
    # scale painting table
    wid<-440 
    ht<-wid*2  
    hm(mat, paste(label,"_DLHM2", sep = ""), wid, ht) # heatmap
  }
  
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
  hmLabel<-paste(mergeLb,"/",mergeLb,"_",gsLabel,sep="") # saving dir included
  # plotting
  matchNa<-exprHeatmap(gsItem, mergeNMdata,hmLabel) # for refining
  lfcHeatmap(gsItem, mergeLFCdata,hmLabel)
  dlfcHeatmap(gsItem, mergeLFCdata,hmLabel)
  
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

