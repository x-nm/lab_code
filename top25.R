# top25 of W-H-PO-OB25
# TOP 25 expressed genes, top 25 upregulated genes, top 25 down regulated genes.
# 
# heatmap functions in figures_func.R
#
# 2017.7.20 by xnm 


wd1<-"D:\\xnm\\work\\17.5.25\\LFC"
wd2<-"D:\\xnm\\work\\17.5.25\\Small Figures/"
saveDir<-"D:\\xnm\\work\\17.5.25\\top25/"

setwd(saveDir)

fdr<-0.1
fc<-0.5

source("D:\\xnm\\work\\CODE\\figures_func.R")
options(stringsAsFactors = FALSE) # must, or ID will be transformed into num...

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                             1. NORM                                                   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
setwd(wd2)
WJHN<-read.table("WJHN_log2_norm_mean_no_aggr.txt",sep = "\t") # norm mean of WJHN
dat_merge<-read.table("WJHNO_log2_norm_mean_aggr.txt",sep = "\t") # get norm_mean of PO

top_W_norm<-top_norm(WJHN,"W","norm_mean",25)
top_J_norm<-top_norm(WJHN,"J","norm_mean")
top_H_norm<-top_norm(WJHN,"H","norm_mean")
top_N_norm<-top_norm(WJHN,"N","norm_mean")
top_P_norm<-top_norm(dat_merge,"P","norm_mean")
top_O_norm<-top_norm(dat_merge,"O","norm_mean")

norm<-c()
sps<-list(top_W_norm,top_J_norm,top_H_norm,top_N_norm,top_O_norm,top_P_norm)
lbs<-c("W","J","H","N","O","P")
for(i in 1:length(sps)){
  dat<-cbind(sps[[i]],rep(lbs[i],nrow(sps[[i]])))
  colnames(dat)<-c("ID","Value","Sample")
  norm<-rbind(norm,dat)
}
norm<-transform(norm,ID = factor(ID, levels = unique(norm$ID)), Sample = factor(Sample, levels = unique(lbs)))
# dat for simple_heatmap: ID, Value, Sample
simple_heatmap(norm,"Top 25 Expressed Genes")
ggsave(paste0(saveDir,"top25_norm_mean_2.png"),width = 4,height = 12)
simple_heatmap(norm,"Top 25 Expressed Genes",vec = F)
ggsave(paste0(saveDir,"top25_norm_mean.png"),width = 12,height = 4)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                         1. LFC: NO aggr                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
setwd(wd1)
W_dat<-read.table("WJA_LFC.txt",sep = "\t")
H_dat<-read.table("HNA_LFC.txt",sep = "\t")
OB25_dat<-read.table("OB25_LFC.txt",sep = "\t")
PO_dat<-read.table("PO_LFC.txt",sep = "\t")

W_DEG<-dat_DEG(W_dat,fdr = fdr,fc = fc)
H_DEG<-dat_DEG(H_dat,fdr = fdr,fc = fc)
OB25_DEG<-dat_DEG(OB25_dat,fdr = fdr,fc = fc)
PO_DEG<-dat_DEG(PO_dat,fdr,fc)

sps<-list(W_DEG,H_DEG,OB25_DEG,PO_DEG)
lbs<-c("WJ","HN","OB25","PO")
up<-c()
down<-c()
for(i in 1:length(lbs)){
  dat1<-top_deg(sps[[i]],up = T,lbs[i],"deg_up")
  dat1<-cbind(dat1,rep(lbs[i],nrow(dat1)))
  colnames(dat1)<-c("ID","Value","Sample")
  up<-rbind(up,dat1)
  
  dat2<-top_deg(sps[[i]],up = F,lbs[i],"deg_down")
  dat2<-cbind(dat2,rep(lbs[i],nrow(dat2)))
  colnames(dat2)<-c("ID","Value","Sample")
  down<-rbind(down,dat2)
}
up<-transform(up, ID=factor(ID, levels = unique(up$ID)), Sample = factor(Sample, levels = unique(lbs)))
down<-transform(down, ID=factor(ID, levels = unique(down$ID)), Sample = factor(Sample, levels = unique(lbs)))

simple_heatmap(up,"Top 25 Up Regulated Genes",vec = T)
ggsave(paste0(saveDir,"top25_up_1.png"),width = 4,height = 12)
simple_heatmap(up,"Top 25 Up Regulated Genes",vec = F)
ggsave(paste0(saveDir,"top25_up_2.png"),width = 12,height = 4)

simple_heatmap(down,"Top 25 Down Regulated Genes",vec = T)
ggsave(paste0(saveDir,"top25_down_1.png"),width = 4,height = 12)
simple_heatmap(down,"Top 25 Down Regulated Genes",vec = F)
ggsave(paste0(saveDir,"top25_down_2.png"),width = 12,height = 4)


# p-val version: ordered by p-val, valued by LFC
up<-c()
down<-c()
for(i in 1:length(lbs)){
  dat1<-top_deg_pval(sps[[i]],up = T,lbs[i],"deg_up_pval", 25)
  dat1<-cbind(dat1,rep(lbs[i],nrow(dat1)))
  colnames(dat1)<-c("ID","Value","Sample")
  up<-rbind(up,dat1)
  
  dat2<-top_deg_pval(sps[[i]],up = F,lbs[i],"deg_down_pval")
  dat2<-cbind(dat2,rep(lbs[i],nrow(dat2)))
  colnames(dat2)<-c("ID","Value","Sample")
  down<-rbind(down,dat2)
}
up<-transform(up, ID=factor(ID, levels = unique(up$ID)), Sample = factor(Sample, levels = unique(lbs)))
down<-transform(down, ID=factor(ID, levels = unique(down$ID)), Sample = factor(Sample, levels = unique(lbs)))

simple_heatmap(up,"Top 25 Up Regulated Genes",vec = T)
ggsave(paste0(saveDir,"top25_up_pval_1.png"),width = 4,height = 12)
simple_heatmap(up,"Top 25 Up Regulated Genes",vec = F)
ggsave(paste0(saveDir,"top25_up_pval_2.png"),width = 12,height = 4)

simple_heatmap(down,"Top 25 Down Regulated Genes",vec = T)
ggsave(paste0(saveDir,"top25_down_pval_1.png"),width = 4,height = 12)
simple_heatmap(down,"Top 25 Down Regulated Genes",vec = F)
ggsave(paste0(saveDir,"top25_down_pval_2.png"),width = 12,height = 4)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                         2. TO GENESET                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

norm_gene<-norm$ID
norm_gene<-norm_gene[!duplicated(norm_gene)]

gsDir<-"D:\\xnm\\work\\17.2.15\\KeyGene\\geneset\\"
write(norm_gene,paste0(gsDir,"top25_norm.gmx"))

up_gene<-up$ID
up_gene<-up_gene[!duplicated(up_gene)]
up_gene<-c("TOP25_UP_LFC",">",up_gene)
write(up_gene,paste0(gsDir,"TOP25_UP_LFC.gmx"))

down_gene<-down$ID
down_gene<-down_gene[!duplicated(down_gene)]
down_gene<-c("TOP25_DOWN_LFC",">",down_gene)
write(down_gene,paste0(gsDir,"TOP25_DOWN_LFC.gmx"))


up_gene<-up$ID
up_gene<-up_gene[!duplicated(up_gene)]
up_gene<-c("TOP25_UP_PVAL",">",up_gene)
write(up_gene,paste0(gsDir,"TOP25_UP_PVAL.gmx"))

down_gene<-down$ID
down_gene<-down_gene[!duplicated(down_gene)]
down_gene<-c("TOP25_DOWN_PVAL",">",down_gene)
write(down_gene,paste0(gsDir,"TOP25_DOWN_PVAL.gmx"))