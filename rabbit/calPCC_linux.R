# 2016.10.27 BY xnm
# INPUT: "normValue_HN.txt"
# OUTPUT: "calPCC_HN.RData"
# usage: nohup Rscript calPCC_linux.R HN >report.txt 2>&1 &

setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\WCGNA\\Data")

rm(list=ls(all=T))


# Args <- commandArgs(trailingOnly = TRUE)
# label <- Args[1]
label <- "HN"



FILE <- paste("normValue_",label,".txt", sep="")

# ##TEST##########################################################
# label2 <-"WJ"
# FILE2 <- paste("normValue_",label2,".txt", sep="")
# 
# dataAll2 <- read.table(FILE2, header = T, row.names = 1)
# dataExpr2 <- dataAll2[,1:spNumExpr] # 前4列为实验组
# dataCtl2 <- dataAll2[,(spNumExpr+1):spNumAll] # 后四列为对照组
# 
# index<-match(row.names(dataCtl),row.names(dataCtl2))
# dataExpr2<-dataExpr2[index,]
# dataExpr2[is.na(dataExpr2)]<-0
# dataCtl2<-dataCtl2[index,]
# dataCtl2[is.na(dataCtl2)]<-0
# 
# dataExpr <- cbind(dataExpr,dataExpr2)
# dataCtl<-cbind(dataCtl,dataCtl2)
# 
# #############################################################

dataAll <- read.table(FILE, header = T, row.names = 1)
spNumAll <- 8
spNumExpr <- 4 # sample number of experimental group
spNumCtl <- 4 # sample number of control
geneNum <- length(rownames(dataAll))

dataExpr <- dataAll[,1:spNumExpr] # 前4列为实验组
dataCtl <- dataAll[,(spNumExpr+1):spNumAll] # 后四列为对照组

# # #TEST
# geneNum <- 100
# dataExpr <- dataExpr[1:100,]
# dataCtl <- dataCtl[1:100,]

delta1 = 0.65
delta2 = 1.5 # try 看boxplot感觉差不多

# 1. Calculate PCC and DE of delta
corValExpr <- round(cor(t(dataExpr), method = "pearson"),4)
corValCtl <- round(cor(t(dataCtl), method = "pearson"),4)

# 2. for delta1
diffCorVal <- abs(corValExpr) -abs(corValCtl)
diffDelta <- abs(diffCorVal) > delta1

# get the (u,v) pairs
uvPairs<-c() #i.e. DE of delta
count<-2
for (i in 1:(geneNum-1)){
  for (j in count:geneNum){
    if(!is.na(diffDelta[i,j])){ # to avoid error due to NA...
      if(diffDelta[i,j]=="TRUE"){
        uvPairs<-rbind(uvPairs,c(i,j))
      }
    }
  }
  count<-count+1
}
colnames(uvPairs)<-c("u","v")

ct1 <-dim(uvPairs)[1]


# 3. for delta2
diffCorVal <- corValExpr - corValCtl
diffDelta <- abs(diffCorVal) > delta2

count<-2
for (i in 1:(geneNum-1)){
  for (j in count:geneNum){
    if(!is.na(diffDelta[i,j])){ # to avoid error due to NA...
      if(diffDelta[i,j]=="TRUE"){
        uvPairs<-rbind(uvPairs,c(i,j))
      }
    }
  }
  count<-count+1
}
ct2<-dim(uvPairs)[1]


# boxplot(abs(diffCorVal))

ct<-rbind(ct1,ct2)
row.names(ct)<-c("ct1","ct2")
colnames(ct)<-"count"


setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DE")
write.table(ct,paste("ct_",label,".txt",sep = ""),quote = F)
save(geneNum, spNumCtl, spNumExpr, dataExpr, dataCtl, uvPairs, file = paste("calPCC_",label, ".RData",sep=""))
