setwd("D:\\Files\\Study\\½»´ó\\¿ÎÌâ\\17.1.4\\data_process\\DataSet")

#############################################################################
# clean WJ HN
#############################################################################
WJ<-read.table("WJ_NORM.txt", header = T, row.names = 1)
# PO<-read.table("PO_NORM.txt", header = T, row.names = 1)
HN<-read.table("HN_NORM.txt", header = T, row.names = 1)

geneId2Sym <- read.csv("gene_id2symbol.txt", header = F, sep = "\t")
indexExclude<-which(geneId2Sym[,2]=="")

WJCL<-cbind(geneId2Sym[,2],WJ[1:41554,1:8])
WJCL<-WJCL[-indexExclude,] # NA excluded, 22682 left
colnames(WJCL)[1]<-"ID"

HNCL<-cbind(geneId2Sym[,2],HN[1:41554,1:8])
HNCL<-HNCL[-indexExclude,] # NA excluded, 22682 left
colnames(HNCL)[1]<-"ID"

# row.names(WJCL)<-geneId2Sym[-indexExclude,2]
# TEMP <- geneId2Sym[-indexExclude,]
# dup <- which(duplicated(TEMP)) # 3340 duplicated geneID
# geneDup <- TEMP[dup,]

write.table(WJCL,"WJ_NORM_CL.txt", quote = F, sep = "\t",row.names = F)
write.table(HNCL,"HN_NORM_CL.txt", quote = F, sep = "\t",row.names = F)

#############################################################################
# logFC FILE
#############################################################################
po<-read.table("DEG_PO.txt",header = T)
wj<-read.table("DEG_WJ.txt",header = T)
hn<-read.table("DEG_HN.txt",header = T)
geneId2Sym <- read.csv("gene_id2symbol.txt", header = F, sep = "\t")

# 17.1.12 VERSION DEG_HN.txt IS WRONG WITH SIGN
# hn[,2]<--hn[,2]
# write.table(hn,"DEG_HN.txt",quote = F, sep="\t")

# wjcl<-cbind(geneId2Sym[,2],wj[1:41554,c(2,6)]) # not the same sequence as geneId2Sym
index<-match(row.names(wj),geneId2Sym[,1])
ID<-geneId2Sym[index,2]
wjcl<-cbind(ID[1:41554],wj[1:41554,c(2,6)])
indexExclude<-which(wjcl[,1]=="")
wjcl<-wjcl[-indexExclude,] # NA excluded, 22682 left
colnames(wjcl)[1]<-"ID"

# hncl<-cbind(geneId2Sym[,2],hn[1:41554,c(2,6)])
index<-match(row.names(hn),geneId2Sym[,1])
ID<-geneId2Sym[index,2]
hncl<-cbind(ID[1:41554],hn[1:41554,c(2,6)])
indexExclude<-which(hncl[,1]=="")
hncl<-hncl[-indexExclude,] # NA excluded, 22682 left
colnames(hncl)[1]<-"ID"

po<-po[,c(2,6)]

write.table(po,"PO_LFC.txt", quote = F, sep = "\t")
write.table(wjcl,"WJ_LFC.txt", quote = F, sep = "\t",row.names = F)
write.table(hncl,"HN_LFC.txt", quote = F, sep = "\t",row.names = F)


#############################################################################
# Merge dataframe for all samples
#############################################################################
WJ<-read.table("WJ_NORM_CL.txt", header = T)
PO<-read.table("PO_NORM.txt", header = T)
HN<-read.table("HN_NORM_CL.txt", header = T)

PO<-cbind(row.names(PO),PO)
colnames(PO)[1]<-"ID"

dup<-duplicated(WJ$ID)
sum(dup)

WJ<-WJ[!dup,]
HN<-HN[!dup,]

TEMP<-merge(WJ,HN, by="ID")  
ALL<-merge(TEMP,PO, by="ID")
write.table(ALL,"WHP_NORM.txt",quote = F, sep = "\t", row.names = F)  

# merge with union set
TEMP<-merge(WJ,HN, by="ID", all = T)  
ALL<-merge(TEMP,PO, by="ID", all = T)
write.table(ALL,"WHPA_NORM.txt",quote = F, sep = "\t", row.names = F)  


# lfc
wj<-read.table("WJ_LFC.txt", header = T)
po<-read.table("PO_LFC.txt",header=T)
hn<-read.table("HN_LFC.txt",header = T)

po<-cbind(row.names(po),po)
colnames(po)[1]<-"ID"

dup<-duplicated(wj$ID)
wj<-wj[!dup,]
dup<-duplicated(hn$ID)
hn<-hn[!dup,]

temp<-merge(wj, hn, by="ID")
all<-merge(temp, po, by="ID")
write.table(all,"WHP_LFC.txt",quote = F, sep = "\t", row.names = F)  

write.table(wj,"WJ_LFC.txt", quote = F, sep = "\t",row.names = F)
write.table(hn,"HN_LFC.txt", quote = F, sep = "\t",row.names = F)

# merge with union set
temp<-merge(wj, hn, by="ID", all=T)
all<-merge(temp, po, by="ID", all = T)
write.table(all,"WHPA_LFC.txt",quote = F, sep = "\t", row.names = F) 
  