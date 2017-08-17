#######################################################################################
# new heatmap for vc
#
# Main operation part
#
# 1. separated data groups
# 2. match gene names (repliecated names allowed)
#       SPECIAL WAY TO SHOW THE REPLICATED ITEMS WITH DIFFERENT TIMES OF REP:
#       FIRST GET THE WHOLE SET: AABBCCD; THEN FILL IT BY SAMPLE, FILL IN BLANK IF NA.
# 3. no clustering of data and sample; show heatmap in the given order
#
# 
# INPUT:
# several group identifiers
# NORM. TXT: (Identifier) ID GROUP1 GROUP2 ...
# LFC.TXT: (Identifier) ID LFC padj
#
# * in LFC FILE, ID NA were excluded.
# 
# 
# 2017.5.26 by xnm
# 
#######################################################################################

setwd("D:\\xnm\\work\\17.5.25\\LFC/")
genesetDir<-"D:\\xnm\\work\\17.2.15\\KeyGene\\geneset/"
gmtDir<-"D:\\xnm\\work\\17.5.25\\genesetSel/" # DEGct file in this dir, too
samples<-readLines("samples.txt") # a list of filename
origData<-read.csv("../ori_whole_data_dup_del.txt",header = T, sep = "\t") # original whole data; generated from dup del, to avoid the influnce of dup:[origWhole(samples)]

source("D:\\xnm\\work\\CODE\\newHeatmap_func.R")


#################################
# ossiSp<-c("WJA","HNA","OB25","PO")
# ossiIndex<-match(ossiSp,samples)

fdr<-0.1
fc<-2

dat_all<-wholeDataSet(fdr,fc,samples, aggr = T) 
# DEGct<-read.csv(paste0(gmtDir,"DEGCt_fdr",fdr,"fc",fc,".txt"),sep = "\t") # DEGct depends on fdr/fc
# prior_prob<-DEGct$degCount/sum(DEGct$degCount)


gmtFile<-"c5.all.v6.0.symbols.gmt"
gmt<-read.csv(paste0(gmtDir,gmtFile),header = F, sep = ",",stringsAsFactors = F)

type<-"gmx"
geneset<-"test" #"CHOLESTEROL_HOMEOSTASIS","BONE(.gmx)"

# save.title<-"test.pdf"
#################################

gs<-getGeneset(geneset,type,"list") # a list of gene
hgt<-HGTpval4samples(origData, samples, gs, fdr, fc) # hypergeometric test

dat_sub<-dataSubset(dat_all, gs, 1) # 1, sort as ID sequence; 2, sort as geneset sequence

sps<-c("WJA","HNA")
sps<-c("WJA","HNA","PO","OB25","OB12","OB8","OB2","CV2","CV8","CV12","CV25")
sps<-c("EMB","PO")
sps<-c("WJA","HNA","OB25","CV25","OB2","CV2","HNL","WJL")
dat_sub<-dataSubset4sps(dat_sub,sps,1)

(plot<-hm(dat_sub))
(plot<-hm(dat_sub,type = "blank", sqr = F))

saveDir<-"D:\\xnm\\work\\17.5.25\\LFC\\output/"
ggsave(paste0(saveDir,geneset,"2.png"),width = 5,height = 13)

# saveDir<-"D:\\xnm\\work\\17.5.25\\top25/"
# ggsave(paste0(saveDir,geneset,"_fc2_n8000.png"),width = 3,height = 13)

(info<-getGeneset(geneset,type,"info"))

# dirtribution plot
degMatchedCt<-sapply(samples,degMatchCt4sp,dat = dat_sub) # acturally, use as.vector(table(dat_sub$sample)) is enough
(dPlot<-distriPlot(prior_prob,degMatchedCt)) # as distriTest and score function only allows ct input




# hmSave(plot,save.title, a = 4, b = 3, save.ht = 3.15)


########################### CKD 
ua<-read.table("D:\\xnm\\work\\DATA\\UA\\UA_LFC.txt",header = T)
ua<-cbind(ua,rep("UA",nrow(ua)))
colnames(ua)[4]<-"sample"
origData_ckd<-rbind(origData,ua)
samples<-c("UA","WJA","HNA","PO","OB25","OB12","OB8","OB2","CV2","CV8","CV12","CV25")

ua<-cbind(ua,rep("",nrow(ua)))
colnames(ua)[5]<-"dup"
dat_all_ckd<-rbind(dat_all,ua)

sps<-c("UA","WJA","HNA")
sps<-c("UA","WJA","HNA","OB25","CV25","PO")
sps<-c("UA","WJA","HNA","PO","OB25","OB12","OB8","OB2","CV2","CV8","CV12","CV25")

geneset<-"KNOCK-OUT-TO-VC" 
gs<-getGeneset(geneset,type,"list") # a list of gene
hgt<-HGTpval4samples(origData_ckd, sps, gs, fdr, fc) # hypergeometric test

dat_sub<-dat_all_ckd[dat_all_ckd$ID%in%gs,]
dat_sub<-dataSubset4sps(dat_sub,sps,1)
(plot<-hm(dat_sub))
(plot<-hm(dat_sub,type = "blank",sqr = F))


# ----------------- FOR SINGAL GENE----------------------------------------
geneset<-"RUNX2" 
gs<-geneset

hgt<-HGTpval4samples(origData, samples, gs, fdr, fc) # hypergeometric test
dat_sub<-dataSubset(dat_all, gs, 1) # 1, sort as ID sequence; 2, sort as geneset sequence

(plot<-hm(dat_sub))
