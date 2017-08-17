# get hgt distribution
# also a new version and new framework of genesetSel
#
# 12min for c5
#
# 2017.7.7 by xnm


setwd("D:\\xnm\\work\\17.5.25\\genesetSel/")
genesetDir<-"D:\\xnm\\work\\17.2.15\\KeyGene\\geneset"
samples<-readLines("../LFC/samples.txt")
origData<-read.csv("../ori_whole_data_dup_del.txt",header = T, sep = "\t") # original whole data; generated from dup del, to avoid the influnce of dup

source("D:\\xnm\\work\\CODE\\genesetSel_func.R")

###################################
fc<-0.5
fdr<-0.1

sp_used<-c("WJA","HNA","OB25","PO")

saveTitle<-"gmx_hgt_8000_fc2.txt"
type<-"gmt"

gmtFile<-"h.all.v6.0.symbols.gmt"
gmt<-read.csv(gmtFile,header = F, sep = ",",stringsAsFactors = F)
###################################

dat_all<-wholeDataSet(fdr,fc,samples, aggr = T)

TIME1<-Sys.time()
res_all<-c()
if(type=="gmt"){
  genesetGroup<-c()
  for(i in 1:nrow(gmt)){
    geneset<-as.vector(strsplit(gmt[i,1],split = "\t")[[1]])
    genesetGroup<-c(genesetGroup,geneset[1])
    geneset<-geneset[-c(1,2)] # geneset(a list) is the only key point for the consequent steps
    ############ main ##########################################################################
    hgt<-HGTpval4samples(origData, sp_used, geneset, fdr, fc) # hypergeometric test
    res<-hgt[,4]
    # res<-hgt[,3]
    ############################################################################################
    res_all<-rbind(res_all,res)
  }
}else if(type=="gmx"){
  # genesetGroup of .gmx in a dir
  genesetGroup<-grep(".gmx$",list.files(genesetDir),value = T)
  for(i in 1:length(genesetGroup)){
    geneset<-gainGeneset(genesetGroup[i]) # geneset(a list) is the only key point for the consequent steps
    ############ main ##########################################################################
    hgt<-HGTpval4samples(origData, sp_used, geneset, fdr, fc) # hypergeometric test
    res<-hgt[,4]
    # res<-hgt[,3]
    ############################################################################################
    res_all<-rbind(res_all,res)
  }
}
colnames(res_all)<-sp_used
row.names(res_all)<-genesetGroup

write.table(res_all,saveTitle,quote = F,sep = "\t")

TIME2<-Sys.time()
print(TIME2-TIME1)


## distribution

par(mfrow=c(2,2))

for (i in 1:4){
  barplot(table(res_all[,i]),main = colnames(res_all)[i])
}


#####################################################################################################
a<-res_all[,4]<0.000001
sum(a)
