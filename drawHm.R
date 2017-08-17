# for input a list and LFC MERGE DATA, and draw a LFC heatmap
# 2017.3.7 by xnm

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#                                DEFINATION                                   #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

setwd("D:\\xnm\\work\\17.2.15\\KeyGene")
# mergeLb<-"WHPA" # input data matrix is WHP_NORM.txt, WHP_LFC.txt

library(pheatmap)
library(gplots)

# make a new dir named mergeLb
if (!dir.exists(mergeLb)){
  dir.create(mergeLb)
}


# test
# gsData<-geneList
# lfcData<-mergeLFCdata


####################################### Functions #############################
# 0.0 find the row index: rows with dist() val == NA
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

# 0.1 return the abs ceiling of floor of a digital numbel, depends on the sign 
dpt<-function(x){
  if(x<0){
    x<-abs(floor(x))
  } 
  else{
    x<-ceiling(x)
  }
  return(x)
}
# 0.2 Draw heatmap; 
# mark=13:heatmap with unclustered rows and rownames; 
# 23: unclustered rows and rownames and display numbers.
# 14, 24: unclusterd rows and cols 
# 04: unclustered, no row name, no values
hm<-function(matrix, mark){
  depth<-max(dpt(min(matrix, na.rm = T)), dpt(max(matrix, na.rm = T)))
  c<-seq(-(depth+0.5),(depth), 0.5)
  len<-length(c)
  col<-bluered(len)
  # plotting
  if(mark==0) pheatmap(matrix, border_color = NA, breaks = c, color = col, show_rownames = F)
  else if(mark==1)   pheatmap(matrix, border_color = NA, breaks = c, color = col) # only row.names
  else if(mark==2) pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T) # TO SHOW THE NA
  else if(mark==13) pheatmap(matrix, border_color = NA, breaks = c, color = col, cluster_rows = F) # only row.names
  else if(mark==23) pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T, cluster_rows = F) # TO SHOW THE NA
  else if(mark==14) pheatmap(matrix, border_color = NA, breaks = c, color = col, cluster_rows = F, cluster_cols = F) # only row.names
  else if(mark==24) pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T, cluster_rows = F, cluster_cols = F) # TO SHOW THE NA
  else if(mark==04) pheatmap(matrix, border_color = NA, breaks = c, color = col, show_rownames = F, cluster_rows = F, cluster_cols = F)
}
# 0.3
ehm<-function(matrix){
  pheatmap(matrix, border_color = NA)
}

# 1. Heatmap of gene expr (log(normValue+1))
exprHeatmap<-function(gsData, exprData){
  # get index
  index<-match(gsData, exprData$ID)
  # na<-sum(is.na(index)) # for refining
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
    ehm(mat) # heatmap
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
    ehm(mat1) # heatmap
    
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
    ehm(mat) # heatmap
  }
  # return(na) # for refining
}

# 2. Heatmap of lfc
# Input: genesetData(gsItem), lfcDataframe, label for naming (saving dir included)
# *gsItem is a vector, with geneID only
# *lfc data: ID,lfc,padj
# Output: heatmap with cluster of genes and cluster of ALL SAMPLES for ONE GENESET
lfcHeatmap<-function(gsData, lfcData, mk){
  # get index
  index<-match(gsData, lfcData$ID)
  na<-sum(is.na(index)) # for refining
  index<-na.omit(index)
  # remake the matrix
  num<-ncol(lfcData)
  del<-seq(3,num,by=2)
  mat<-lfcData[index,-del] # exclude the padj col
  row.names(mat)<-mat$ID
  mat<-mat[,-1]
  mat<-as.matrix(mat)
  
  # delete all NA rows
  mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
  
  # draw the heatmap
  if(all(!is.na(dist(mat)))){
    hm(mat,mk) # heatmap
    print("lhm0")
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
    hm(mat1,mk) # heatmap
    print("lhm1")
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
    hm(mat,mk) # heatmap
    print("lhm2")
  }
  return(na) # for refining
  
}

# 2.2 Heatmap of dlfc: with non-DEG colored as gray
# return: DEG DATA MATRIX
dlfcHeatmap<-function(gsData, lfcData, mk){
  # get index
  index<-na.omit(match(gsData, lfcData$ID))
  # remake the matrix
  num<-ncol(lfcData)
  del<-seq(3,num,by=2)
  mat<-lfcData[index,-del] # exclude the padj col
  
  del<-seq(2,num-1,by=2)
  padj<-lfcData[index,-del]
  mat[padj>0.05]<-NA
  # padj[is.na(padj)]<-0
  row.names(mat)<-mat$ID
  mat<-mat[,-1]
  mat<-as.matrix(mat)

  # delete all NA rows
  mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
  res<-mat # DEG DATA
  # draw the heatmap
  if(all(!is.na(dist(mat)))){
    hm(mat,mk) # heatmap
    print("dlhm0")
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
    hm(mat1,mk) # heatmap
    print("dlhm1")
    
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
    # scale painting table
    hm(mat,mk) # heatmap
    print("dlhm2")
  }
  return(res)
}
###############################################################################



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#                                   Main                                      #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

# data
mergeLb<-"WHPOCUa" 
# OBCVa, WHPOCa,WHmOCa,WHm,WHPOCUa
# glLabel<-"MK"
mergeLFCdata<-read.table(paste(mergeLb,"_LFC.txt",sep = ""),header = T)
# mergeNMdata<-read.table(paste(mergeLb,"_NORM.txt", sep = ""),header = T)

# WH-OB:
# index<-c(1,2,3,4,5,8,9,10,11,12,13,14,15)
# mergeLFCdata<-mergeLFCdata[,index]

# geneList
glLabel<-"KG"
# glLabel<-"GO_ECM"
geneList<-readLines(paste("geneset/",glLabel,".gmx",sep = ""))
geneList<-geneList[-c(1,2)]
# DIY
# glLabel<-"test"
# geneList<-c("HDAC3","GPNMB","CRP","P21","AXIN2","SFRP1","PMEL")

mna<-c()
# plotting
# exprHeatmap(geneList, mergeNMdata) # for refining

# 0:hm without num and row.names;1:hm with row.names only ;2:hm with num and row.names; 
matchNa<-lfcHeatmap(geneList, mergeLFCdata,1)
data<-dlfcHeatmap(geneList, mergeLFCdata,1)

# match NA
mna<-rbind(mna,c(glLabel,matchNa,length(geneList)))

# save combined DEG List
write.table(row.names(data),
            paste("geneset/",mergeLb,"_",glLabel,"_DEG.gmx",sep = ""),
            quote = F,sep = "\t",row.names = F, col.names = F)



# For saving:
hmLabel<-paste(mergeLb,"/",mergeLb,"_",glLabel,sep="") # saving dir included
# lfcHeatmap.save(geneList, mergeLFCdata, hmLabel)
dlfcHeatmap.save(geneList, mergeLFCdata, hmLabel)



####################################### Functions for save #############################

hm.save<-function(matrix, title, width, height){
  depth<-max(dpt(min(matrix, na.rm = T)), dpt(max(matrix, na.rm = T)))
  c<-seq(-(depth+0.5),(depth), 0.5)
  len<-length(c)
  col<-bluered(len)
  
  # bmp(paste(title,".bmp",sep=""), width = width, height = height)
  # # plotting
  # pheatmap(matrix, border_color = NA, breaks = c, color = col)
  # dev.off()
  
  bmp(paste(title,"_num.bmp",sep=""), width = width, height = height)
  # plotting
  pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T) # TO SHOW THE NA
  dev.off()
}
ehm.save<-function(matrix, title, width, height){
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
exprHeatmap.save<-function(gsData, exprData, label){
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
    ehm.save(mat, paste(label,"_EHM", sep = ""), wid, ht) # heatmap
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
    ehm.save(mat1, paste(label,"_EHM1", sep = ""), wid, ht) # heatmap
    
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
    ehm.save(mat, paste(label,"_EHM2", sep = ""), wid, ht) # heatmap
  }
  
  return(na) # for refining
}

# 2. Heatmap of lfc
# Input: genesetData(gsItem), lfcDataframe, label for naming (saving dir included)
# *gsItem is a vector, with geneID only
# *lfc data: ID,lfc,padj
# Output: heatmap with cluster of genes and cluster of ALL SAMPLES for ONE GENESET
lfcHeatmap.save<-function(gsData, lfcData, label){
  # get index
  index<-na.omit(match(gsData, lfcData$ID))
  # remake the matrix
  num<-ncol(lfcData)
  del<-seq(3,num,by=2)
  mat<-lfcData[index,-del] # exclude the padj col
  row.names(mat)<-mat$ID
  mat<-mat[,-1]
  mat<-as.matrix(mat)
  
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
    hm.save(mat, paste(label,"_LHM", sep = ""), wid, ht) # heatmap
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
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
    hm.save(mat1, paste(label,"_LHM1", sep = ""), wid, ht) # heatmap
    
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
    hm.save(mat, paste(label,"_LHM2", sep = ""), wid, ht) # heatmap
  }
  
}

# 2.2 Heatmap of dlfc: with non-DEG colored as gray
dlfcHeatmap.save<-function(gsData, lfcData, label){
  # get index
  index<-na.omit(match(gsData, lfcData$ID))
  # remake the matrix
  num<-ncol(lfcData)
  del<-seq(3,num,by=2)
  mat<-lfcData[index,-del] # exclude the padj col
  
  del<-seq(2,num-1,by=2)
  padj<-lfcData[index,-del]
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
    hm.save(mat, paste(label,"_DLHM", sep = ""), wid, ht) # heatmap
  }
  else{
    # 1. plot with NA set as 0
    mat1<-mat
    mat1[is.na(mat1)]<-0
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
    hm.save(mat1, paste(label,"_DLHM1", sep = ""), wid, ht) # heatmap
    
    # 2. plot with rows deleted
    mat<-matDelRows(mat)
    # scale painting table
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
    hm.save(mat, paste(label,"_DLHM2", sep = ""), wid, ht) # heatmap
  }
  
}




########################################### for 3 cols of LFC #######################
# # 2. Heatmap of lfc
# # Input: genesetData(gsItem), lfcDataframe, label for naming (saving dir included)
# # *gsItem is a vector, with geneID only
# # *lfc data: ID,lfc,padj
# # Output: heatmap with cluster of genes and cluster of ALL SAMPLES for ONE GENESET
# lfcHeatmap.save<-function(gsData, lfcData, label){
#   # get index
#   index<-na.omit(match(gsData, lfcData$ID))
#   # remake the matrix
#   num<-ncol(lfcData)
#   del<-seq(3,num,by=2)
#   mat<-lfcData[index,-del] # exclude the padj col
#   row.names(mat)<-mat$ID
#   mat<-mat[,-1]
#   mat<-as.matrix(mat)
#   
#   # delete all NA rows
#   mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
#   
#   # draw the heatmap
#   if(all(!is.na(dist(mat)))){
#     # scale painting table
#     wid<-440 
#     ht<-wid*2  
#     hm.save(mat, paste(label,"_LHM", sep = ""), wid, ht) # heatmap
#   }
#   else{
#     # 1. plot with NA set as 0
#     mat1<-mat
#     mat1[is.na(mat1)]<-0
#     # scale painting table
#     wid<-440 
#     ht<-wid*2  
#     hm.save(mat1, paste(label,"_LHM1", sep = ""), wid, ht) # heatmap
#     
#     # 2. plot with rows deleted
#     mat<-matDelRows(mat)
#     # scale painting table
#     wid<-440 
#     ht<-wid*2  
#     hm.save(mat, paste(label,"_LHM2", sep = ""), wid, ht) # heatmap
#   }
#   
# }
# 
# # 2.2 Heatmap of dlfc: with non-DEG colored as gray
# dlfcHeatmap.save<-function(gsData, lfcData, label){
#   # get index
#   index<-na.omit(match(gsData, lfcData$ID))
#   # remake the matrix
#   num<-ncol(lfcData)
#   del<-seq(3,num,by=2)
#   mat<-lfcData[index,-del] # exclude the padj col
#   
#   del<-seq(2,num-1,by=2)
#   padj<-lfcData[index,-del]
#   mat[padj>0.05]<-NA
#   # padj[is.na(padj)]<-0
#   row.names(mat)<-mat$ID
#   mat<-mat[,-1]
#   mat<-as.matrix(mat)
#   
#   # delete all NA rows
#   mat<-mat[apply(mat,1,function(x) !all(is.na(x))),]
#   
#   # draw the heatmap
#   if(all(!is.na(dist(mat)))){
#     # scale painting table
#     wid<-440 
#     ht<-wid*2  
#     hm.save(mat, paste(label,"_DLHM", sep = ""), wid, ht) # heatmap
#   }
#   else{
#     # 1. plot with NA set as 0
#     mat1<-mat
#     mat1[is.na(mat1)]<-0
#     # scale painting table
#     wid<-440 
#     ht<-wid*2  
#     hm.save(mat1, paste(label,"_DLHM1", sep = ""), wid, ht) # heatmap
#     
#     # 2. plot with rows deleted
#     mat<-matDelRows(mat)
#     # scale painting table
#     wid<-440 
#     ht<-wid*2  
#     hm.save(mat, paste(label,"_DLHM2", sep = ""), wid, ht) # heatmap
#   }
#   
# }
