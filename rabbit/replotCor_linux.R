# REPLOT COR
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

# vari%
naNum1 <- sum(is.na(indexInDEG1))
varNum1 <- geneNum - naNum1
naNum2 <- sum(is.na(indexInDEG2))
varNum2 <- geneNum - naNum2
variPct1 <- varNum1/geneNum
variPct2 <- varNum2/geneNum
variNum <- rbind(c(varNum1,variPct1),c(varNum2,variPct2))
row.names(variNum) <- c("WJ","HN")
colnames(variNum) <- c("variNum","percent") 

# combine matrix
lgFCMat <- cbind(DEG1[indexInDEG1,1],DEG2[indexInDEG2,1])
colnames(lgFCMat) <- c("WJ","HN")
row.names(lgFCMat) <- GENELIST
lgFCMat[is.na(lgFCMat)]<-0


# plot.1 corr plot of logFC ##########################################################
title <- paste("corr_",lable,"_logFC", sep = "")
x <- lgFCMat[,1]
y <- lgFCMat[,2]
fit <- lm(y ~ x)
R.sqr = summary(fit)$adj.r.squared 

# use ggplot2
library(ggplot2)
xLoc <- min(x)+0.85*(max(x)-min(x))
yLoc <- min(y)+0.1*(max(y)-min(y))
ggplot(as.data.frame(lgFCMat), aes(x = x, y = y))+geom_point()+geom_smooth(method=lm,color="red")+
  labs(x="WHHL log(FC)",
       y="NZW-HC log(FC)",
       title=paste("Correlation of Expression Changes of HW in ",lable, sep = ""))+
  annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))
ggsave(paste(title,".png",sep=""),width = 5, height = 5)
