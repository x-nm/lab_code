#2017.3.8 by xnm

# WHm + OCa
w<-read.table("WJ_LFC.txt", header = T)
h<-read.table("HN_LFC.txt",header = T)
wh<-merge(w, h, by = "ID",all = T)
colnames(wh)<-c("ID","WJ","W.padj","HN","H.padj")
write.table(wh,"WHm_LFC.txt", quote = F, sep = "\t", row.names = F)
colnames(wh)<-c("ID","WJ","W.padj","HN","H.padj")
oc<-read.table("OBCVa_LFC.txt",header = T)
M<-merge(wh,oc,by = "ID", all = T)
write.table(M,"WHmOCa_LFC.txt", quote = F, sep = "\t", row.names = F)



#####################################################################################
# GET WJa_LFC, HNa_LFC 还要改代码
#####################################################################################
setwd("D:\\xnm\\work\\DATA\\WH")
# # count item numbers under different padj thresholds
# padjCt<-function(data){
#   return(c(nrow(subset(data, data$pdaj<0.01)),
#            nrow(subset(data, data$pdaj<0.05)),
#            nrow(subset(data, data$pdaj<0.1))))
# }
# geneId2Sym <- read.csv("gene_id2symbol.txt", header = F, sep = "\t")
# indexExclude<-which(geneId2Sym[,2]=="")
# # TO GET THE logFC FILE FROM DEG FILE
# label<-"WJ"
# lfc<-read.table(paste("DEG_",label,".txt",sep = ""),header = T)
# lfc<-cbind(geneId2Sym[,2],lfc[1:41554,c(2,6)])
# lfc<-lfc[-indexExclude,] # NA excluded, 22682 left
# 
# colnames(lfc)<-c("ID","LFC","padj")
# count<-c()
# count<-rbind(count,c("state","<0.01","<0.05","<0.1"))
# count<-rbind(count,c("Before",padjCt(lfc)))
# 
# # clean: deal with duplicated ID
# dup<-which(duplicated(lfc$ID))
# dupItems<-lfc[dup,]
# # lfc_m: just minus the duplicated items 
# lfc_m<-lfc[-dup,]
# # lfc_a: add up duplicted items' value
# lfc_a<-lfc_m
# mch<-match(dupItems$ID,lfc_m$ID)
# for(i in 1:length(dup)){
#   a<-lfc_m[mch[i],]$LFC #用lfc_a则反复相加，但可能偏差较大，再考虑到padj，先用只加前两项？但这样是后两项啊...
#   b<-dupItems[i,]$LFC
#   if(!is.na(a)){
#     if(!is.na(b)) lfc_a[mch[i],]$LFC<-a+b
#   }
#   else{
#     if(!is.na(b)) lfc_a[mch[i],]$LFC<-b
#   }
# }
# # SHOULD CARE ABOUT PADJ?
# count<-rbind(count,c("After",padjCt(lfc_a)))
# 
# write.table(lfc_m,paste(label,"m_LFC.txt", sep = ""),quote = F, row.names = F, sep = "\t")
# write.table(lfc_a,paste(label,"a_LFC.txt", sep = ""),quote = F, row.names = F, sep = "\t")
# write.table(count,paste(label,"_count.txt", sep = ""),quote = F, row.names = F, sep = "\t", col.names = F)
# 
# sum(is.na(dupItems$LFC))
# sum(is.na(lfc_m[mch,]$LFC))
# sum(is.na(lfc_a[mch,]$LFC))
# all(lfc_m[mch,]$ID==dupItems$ID)
# temp<-lfc_a

#####################################################################################


# U<-read.table("UA_LFC.txt",header = T)
# colnames(U)<-c("ID","UA","padj.u")
# 
# A<-read.table("WHPOCa_LFC.txt", header = T)
# M<-merge(A,U,by = "ID", all = T)
# write.table(M,"WHPOCUa_LFC.txt", quote = F, sep = "\t", row.names = F)


# setwd("D:\\xnm\\work\\17.2.15\\KeyGene")
# whp<-read.table("WHPA_LFC.txt",header = T)
# oc<-read.table("OBCVa_LFC.txt",header = T)
# M<-merge(whp,oc,by = "ID", all = T)
# write.table(M,"WHPOCa_LFC.txt", quote = F, sep = "\t", row.names = F)

# whp<-read.table("WHPA_LFC.txt",header = T)
# oc<-read.table("OBCVm_LFC.txt",header = T)
# M<-merge(whp,oc,by = "ID", all = T)
# write.table(M,"WHPOCm_LFC.txt", quote = F, sep = "\t", row.names = F)

# # OB<-read.table("OBm_LFC.txt",header = T)
# # CV<-read.table("CVm_LFC.txt",header = T)
# # m<-merge(OB, CV, by="ID")
# # write.table(m,"OBCVm_LFC.txt", quote = F, sep = "\t", row.names = F)

# OB<-read.table("OBa_LFC.txt",header = T)
# CV<-read.table("CVa_LFC.txt",header = T)
# m<-merge(OB, CV, by="ID")
# write.table(m,"OBCVa_LFC.txt", quote = F, sep = "\t", row.names = F)


#########################################################################################
# OB2<-read.table("OB2a_LFC.txt", header = T)
# OB8<-read.table("OB8a_LFC.txt", header = T)
# OB12<-read.table("OB12a_LFC.txt", header = T)
# OB25<-read.table("OB25a_LFC.txt", header = T)
# M1<-merge(OB2,OB8, by="ID")
# M2<-merge(OB12,OB25, by="ID")
# M<-merge(M1,M2,by="ID")
# colnames(M)<-c("ID","OB2","padj.2","OB8","padj.8","OB12","padj.12","OB25","padj.25")
# write.table(M,"OBa_LFC.txt", quote = F, sep = "\t", row.names = F)

# OB2<-read.table("OB2m_LFC.txt", header = T)
# OB8<-read.table("OB8m_LFC.txt", header = T)
# OB12<-read.table("OB12m_LFC.txt", header = T)
# OB25<-read.table("OB25m_LFC.txt", header = T)
# M1<-merge(OB2,OB8, by="ID")
# M2<-merge(OB12,OB25, by="ID")
# M<-merge(M1,M2,by="ID")
# colnames(M)<-c("ID","OB2","padj.2","OB8","padj.8","OB12","padj.12","OB25","padj.25")
# write.table(M,"OBm_LFC.txt", quote = F, sep = "\t", row.names = F)

# CV2<-read.table("CV2a_LFC.txt", header = T)
# CV8<-read.table("CV8a_LFC.txt", header = T)
# CV12<-read.table("CV12a_LFC.txt", header = T)
# CV25<-read.table("CV25a_LFC.txt", header = T)
# M1<-merge(CV2,CV8, by="ID")
# M2<-merge(CV12,CV25, by="ID")
# M<-merge(M1,M2,by="ID")
# colnames(M)<-c("ID","CV2","padj.2","CV8","padj.8","CV12","padj.12","CV25","padj.25")
# write.table(M,"CVa_LFC.txt", quote = F, sep = "\t", row.names = F)

# CV2<-read.table("CV2m_LFC.txt", header = T)
# CV8<-read.table("CV8m_LFC.txt", header = T)
# CV12<-read.table("CV12m_LFC.txt", header = T)
# CV25<-read.table("CV25m_LFC.txt", header = T)
# M1<-merge(CV2,CV8, by="ID")
# M2<-merge(CV12,CV25, by="ID")
# M<-merge(M1,M2,by="ID")
# colnames(M)<-c("ID","CV2","padj.2","CV8","padj.8","CV12","padj.12","CV25","padj.25")
# write.table(M,"CVm_LFC.txt", quote = F, sep = "\t", row.names = F)
