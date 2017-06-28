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
origData<-read.csv("../fdr1fc0.txt",header = T, sep = "\t") # original whole data

source("D:\\xnm\\work\\CODE\\newHeatmap_func.R")


#################################
ossiSp<-c("WJA","HNA","OB25","PO")
ossiIndex<-match(ossiSp,samples)

fdr<-0.1
fc<-0.5

dat_all<-wholeDataSet(fdr,fc,samples, aggr = T) 
DEGct<-read.csv(paste0(gmtDir,"DEGCt_fdr",fdr,"fc",fc,".txt"),sep = "\t") # DEGct depends on fdr/fc
prior_prob<-DEGct$degCount/sum(DEGct$degCount)


gmtFile<-"c5.all.v6.0.symbols.gmt"
gmt<-read.csv(paste0(gmtDir,gmtFile),header = F, sep = ",",stringsAsFactors = F)

type<-"gmx"
geneset<-"MARKERS-BONE-zeugopod-mesenchymal-condensate-cells" #"CHOLESTEROL_HOMEOSTASIS","BONE(.gmx)"

# save.title<-"test.pdf"
#################################

gs<-getGeneset(geneset,type,"list") # a list of gene
hgt<-HGTpval4samples(origData, samples, gs, fdr, fc) # hypergeometric test

dat_sub<-dataSubset(dat_all, gs, 1) # 1, sort as ID sequence; 2, sort as geneset sequence
(plot<-hm(dat_sub))
(plot<-hm(dat_sub,type = "blank", sqr = F))

# dirtribution plot
degMatchedCt<-sapply(samples,degMatchCt4sp,dat = dat_sub) # acturally, use as.vector(table(dat_sub$sample)) is enough
(dPlot<-distriPlot(prior_prob,degMatchedCt)) # as distriTest and score function only allows ct input



(info<-getGeneset(geneset,type,"info"))

# hmSave(plot,save.title, a = 4, b = 3, save.ht = 3.15)

