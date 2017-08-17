# modified and corrected "otherTissuDataPrepare.R"
# RERUN for the correct LFC, AND also get the LFC file for HM
# 2017.07.18 by xnm


setwd("D:\\xnm\\work\\Rabbit_RNASeq")
saveDir<-"D:\\xnm\\work\\DATA\\WH\\deseq2_output/"
library("DESeq2")
options(stringsAsFactors = FALSE)

raw<-read.table("raw_count.txt",header = T)
colnames(raw)[1]<-"67" # colname of "67" will be turned into "X67"
sps<-read.table("samples.txt",header = F)
geneId2Sym <- read.csv("gene_id2symbol.txt", header = F, sep = "\t")

# Functions:
# 0.1 From symbol abbr to long symbol:
sybTrans <- function(symbol){
  if (symbol=="W" ) return("WHHL")
  if (symbol=="J" ) return("JW")
  if (symbol=="N" ) return("NZW")
  if (symbol=="H" ) return("NZW-HC")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Queries
queries<-c("HN_Aorta","HN_Heart","HN_Liver","HN_Kidney","HN_Heart_coronary",
           "WJ_Kidney","WJ_Aorta","WJ_Heart","WJ_Liver","NZW_embryo_VS_aorta")

# Main
for(i in 1:length(queries)){
  query<-queries[i]
  
  # 1. get labels/index of query sample data group
  if(query=="NZW_embryo_VS_aorta"){
    querySp1Sps<-subset(sps, sps[2]=="NZW-embryo")[,1]
    querySp2Sps<-subset(sps, sps[2]=="NZW-HC-Aorta")[,1]
    saveTitle<-query
    saveTitle2<-"EMB"
  }else if(query=="HN_Heart_coronary"){
    querySp1Sps<-subset(sps, sps[2]=="NZW-HC-Heart/coronary")[,1]
    querySp2Sps<-subset(sps, sps[2]=="NZW-Heart/coronary")[,1]
    saveTitle<-query
    saveTitle2<-"HNHC"
  }else{
    sp1<-substr(query,1,1) # CASE
    sp2<-substr(query,2,2) # CTL
    tissue<-strsplit(query,"_")[[1]][2]
    saveTitle<-paste(sp1,sp2,"_",tissue, sep="") # for save the orig LFC data
    saveTitle2<-paste0(sp1,sp2,substr(tissue,1,1))
    sp1<-sybTrans(sp1)
    sp2<-sybTrans(sp2)
    querySp1Sps<-subset(sps, sps[2]==paste0(sp1,"-",tissue))[,1]
    querySp2Sps<-subset(sps, sps[2]==paste0(sp2,"-",tissue))[,1]
  }
  
  # 2. data prepare
  sp1Data <- raw[,match(querySp1Sps, colnames(raw))]
  sp2Data <- raw[,match(querySp2Sps, colnames(raw))]
  countData <- cbind(sp1Data, sp2Data)
  num <- ncol(sp1Data)
  sp1Col <- cbind(names(sp1Data),c(rep(1,num)))
  sp2Col <- cbind(names(sp2Data),c(rep(0,num)))
  colData <- rbind(sp1Col, sp2Col) # whether the same format as former?
  colnames(colData) <- c("id", "condition")
  
  # 3. DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  # 4. get orig LFC
  # ID  GENE  baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
  resOrdered <- res[order(res$padj),]
  resOrdered <- as.data.frame(resOrdered)
  index<-match(row.names(resOrdered),geneId2Sym[,1]) # 5 NAs in this step, and the NAs are not in the last 5 positions,
  ID<-geneId2Sym[index,2]
  LFC<-cbind(ID,resOrdered) 
  LFC<-LFC[!is.na(LFC$ID),]
  
  indexExclude<-which(LFC[,1]=="")
  LFC<-LFC[-indexExclude,] # ""excluded
  colnames(LFC)[1]<-"ID"
  write.table(LFC,file=paste0(saveDir,saveTitle, "_LFC.txt"), sep = "\t", quote = F)
  
  # 5. get LFC for heatmap input
  dat<-LFC[,c(1,3,7)]
  colnames(dat)<-c("ID","LFC","padj")
  write.table(dat, paste0(saveDir,"LFC_new/",saveTitle2,"_LFC.txt"), sep = "\t", quote = F)
  
}

