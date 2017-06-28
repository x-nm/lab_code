#######################################################################################
# other tissue, raw -> DEG -> LFC(WITH PVAL)+NORM,annotated.
#
# Input:
#    gene_id2symbol.txt
#    raw_count.txt
#    samples.txt
#    query: WJ_Kidney
# Output:
#    group_tissue_LFC.txt
#    group_tissue_NORM.txt
#
# 2017.5.10 by xnm
# 
#######################################################################################


setwd("D:\\xnm\\work\\Rabbit_RNASeq")
library("DESeq2")
options(stringsAsFactors = FALSE)

raw<-read.table("raw_count.txt",header = T)
colnames(raw)[1]<-"67" # colname of "67" will be turned into "X67"
sps<-read.table("samples.txt",header = F)

# Functions:
# 0.1 From symbol abbr to long symbol:
sybTrans <- function(symbol){
  if (symbol=="W" ) return("WHHL")
  if (symbol=="J" ) return("JW")
  if (symbol=="N" ) return("NZW")
  if (symbol=="H" ) return("NZW-HC")
}


#######################################################################################
# 1. Data load and choose, selectively
#######################################################################################

# DEFINE THE QUERIES IN THIS OPERATION


#######################################################################################
query<-"HN_Aorta" # Heart/coronary(HN),WJ_Heart, WJ_Liver, HN_Kidney
#######################################################################################



sp1<-substr(query,1,1) # CASE
sp2<-substr(query,2,2) # CTL
tissue<-strsplit(query,"_")[[1]][2]
# If input by define sp1, sp2, tissue directly:
query<-paste(sp1,sp2,"_",tissue, sep="")
sp1<-sybTrans(sp1)
sp2<-sybTrans(sp2)
# FOR Heart/coronary, in case not trouble saving results:
# query<-"HN_Heart_coronary" 

querySp1Sps<-subset(sps, sps[2]==paste0(sp1,"-",tissue))[,1]
# the sps get from samples file are in different sequence with those in raw data.
querySp2Sps<-subset(sps, sps[2]==paste0(sp2,"-",tissue))[,1]


################################
# # FOR NZW-embryo VS NZW-HC-Aorta
# query<-"NZW_embryo_VS_aorta"
# querySp1Sps<-subset(sps, sps[2]=="NZW-embryo")[,1]
# querySp2Sps<-subset(sps, sps[2]=="NZW-HC-Aorta")[,1]
################################

sp1Data <- raw[,match(querySp1Sps, colnames(raw))]
sp2Data <- raw[,match(querySp2Sps, colnames(raw))]
countData <- cbind(sp1Data, sp2Data)
num <- ncol(sp1Data)
sp1Col <- cbind(names(sp1Data),c(rep(1,num)))
sp2Col <- cbind(names(sp2Data),c(rep(0,num)))
colData <- rbind(sp1Col, sp2Col) # whether the same format as former?
colnames(colData) <- c("id", "condition")

#######################################################################################

# the rest parts of code
source("D:\\xnm\\work\\CODE\\otherTissueDataPrepare_script2.R")

