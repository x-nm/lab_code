# expr heatmap
# 2017.07.07 by xnm

setwd("D:\\xnm\\work\\17.5.25\\LFC/")


source("D:\\xnm\\work\\CODE\\figures_1_func.R")
options(stringsAsFactors = FALSE) # must, or ID will be transformed into num...

################################
fdr<-0.1
fc<-0.5
samples<-readLines("samples.txt") # a list of filename
dat_all<-wholeDataSet(fdr,fc,samples, aggr = T)
################################

##############################################################################################
# 0.1 a general overview for the express profile of all the data used.
# to show the differences.（只作DEG并集的）

# library(reshape2)

# 0.1.1 LFC of all data: ggplot2
gs<-dat_all$ID[!duplicated(dat_all$ID)] # the union of all DEGs, 18191 items for fdr0.1 fc0.5
dat_sub<-dataSubset(dat_all, gs, 1) # 1, sort as ID sequence; 2, sort as geneset sequence
geneset<-"LFC of all data" # title of heatmap
(plot<-hm(dat_sub,type = "blank", sqr = F, hgtStar = F)) # not pretty (too many dup *, and not watchable)

# 0.1.2 LFC of all data: pheatmap
dat_all<-transform(dat_all, sample=factor(sample, levels = unique(samples)))
c<-acast(dat_all,ID~sample, value.var = "LFC")
hm(c,04) # pheatmap, cluster is not available for so many NAs
#
sps<-c("WJA","HNA","OB25","PO","HNL")
dat<-dat_all[dat_all$sample%in%sps,]
dat<-transform(dat, sample=factor(sample, levels = unique(sps)))
c2<-acast(dat,ID~sample, value.var = "LFC")

hm(c2,04)

# 0.1.3 expr of all data
# we got a version of WHP, in 2.15
# only WH+P have norm value


normdir<-"D:\\xnm\\work\\DATA\\NORM"
normfiles<-list.files(normdir)
dat_all<-read.table(paste0(normdir,"\\",normfiles[1]),header=T)
dat_all<-aggrDup(dat_all)
for (i in 2:length(normfiles)){
  f<-normfiles[i]
  dat<-read.table(paste0(normdir,"\\",f),header=T)
  # IT SEEMS MERGE IS NOT AFFECTED BY DUP,BUT WILL MULTIPLY MUCH FOLDS...
  dat<-aggrDup(dat)
  dat_all<-merge(dat_all,dat,by="ID",all=T)
}

hm(dat[,-1],04)

# even not able to draw one with 8 col
# Error: cannot allocate vector of size 470.4 Mb
