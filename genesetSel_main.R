# genesetSel MAIN part
# for a collection of genesets, select genesets of interest
#
# IMPORTANT: samples are collected in the order of samples.txt in LFC dir
#
# note: 
#   "geneset" means a list of genes here, not a filename
#   it can be gained from a gainGeneset() function, from the filename
#
# 2017.6.7 by xnm


setwd("D:\\xnm\\work\\17.5.25\\genesetSel/")
genesetDir<-"D:\\xnm\\work\\17.2.15\\KeyGene\\geneset"
samples<-readLines("../LFC/samples.txt")

source("D:\\xnm\\work\\CODE\\genesetSel_func.R")

###################################
fc<-0.5
fdr<-0.1
ossiSp<-c("WJA","HNA","OB25","PO")
saveTitle<-"h.txt"
type<-"gmx"

gmtFile<-"h.all.v6.0.symbols.gmt"
gmt<-read.csv(gmtFile,header = F, sep = ",",stringsAsFactors = F)
###################################

DEGct<-read.csv(paste0("DEGCt_fdr",fdr,"fc",fc,".txt"),sep = "\t") # DEGct depends on fdr/fc
prior<-DEGct$degCount/sum(DEGct$degCount)
dat_all<-wholeDataSet(fdr,fc,samples)
ossiIndex<-match(ossiSp,samples)

# get the result
if(type=="gmt"){
  # .gmt, get the res_all for a .gmt
  source("D:\\xnm\\work\\CODE\\genesetSel_res_gmt.R") 
}else if(type=="gmx"){
  # genesetGroup of .gmx in a dir
  genesetGroup<-grep(".gmx$",list.files(genesetDir),value = T)
  source("D:\\xnm\\work\\CODE\\genesetSel_res_gmx.R") 
}

write.table(res_all,saveTitle,quote = F,sep = "\t")






#### Followed analysis #########################################################
# res<-read.csv(saveTitle,header = T, sep = "\t")
a<-grep("BONE",row.names(c5),value = T) # grep "OSTEO"
a<-c5[a,]
