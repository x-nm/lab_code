# reheatmap of lgFC withou na
# 2016.11.07

options(stringsAsFactors = FALSE);

Args <- commandArgs(trailingOnly = TRUE)
lable <- Args[1]
FILE <- Args[2]

#####################################################################################


# GENE LIST (用read.table读取所以用的时候要像矩阵一样用才行...)
GENELIST <- read.table(FILE, header = F)
names(GENELIST)<-"GENE"
GENELIST <- GENELIST$GENE #变成向量
geneNum <- length(GENELIST)

# input data
geneId2Sym <- read.csv("gene_id2symbol.txt", header = F, sep = "\t")
names(geneId2Sym)<-c("ID","GeneSymbol")
DEG1 <- read.csv("DEG_WJ.txt",header = T,sep = "\t")# WJ as 1, HN as 2
DEG1 <- DEG1[,c(2,7)]
DEG2 <- read.csv("DEG_HN.txt",header = T,sep = "\t")
DEG2 <- DEG2[,c(2,7)]



###Functions#########################################################################
# get index
indexInDEG1 <- match(GENELIST,DEG1$sy) 
indexInDEG2 <- match(GENELIST,DEG2$sy)
index <- match(GENELIST,geneId2Sym$GeneSymbol)
geneID <- geneId2Sym[index,1]

# combine matrix
lgFCMat <- cbind(DEG1[indexInDEG1,1],DEG2[indexInDEG2,1])
colnames(lgFCMat) <- c("WJ","HN")
row.names(lgFCMat) <- GENELIST
lgFCMat[is.na(lgFCMat)]<-0
index<-rowSums(lgFCMat)!=0
# scaleNum<-length(index)/53

library(pheatmap)

title <- paste("HM_",lable,"_logFC_nona", sep = "")
bmp(paste(title,".bmp",sep=""), width = 500, height = 960)
pheatmap(lgFCMat[index,],fontsize = 11)
dev.off()
