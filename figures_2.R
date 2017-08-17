# log(norm+1) corr plot
# 2017.7.11 by xnm


setwd("D:\\xnm\\work\\17.5.25\\Small Figures/")


source("D:\\xnm\\work\\CODE\\figures_func.R")
options(stringsAsFactors = FALSE) # must, or ID will be transformed into num...

#################################
# fdr<-0.1
# fc<-0.5
# samples<-readLines("samples.txt") # a list of filename
# dat_all<-wholeDataSet(fdr,fc,samples, aggr = T) 
#################################

# 0.2 expression levels dot plot of W-J, H-N, W-H, W-OB25, H-OB25. W-PO?, H-PO?
# point out genes with high expression values in the diagonal. calculate pearson's R , R^2
# log(x+1) transformation??¡¾¡¿
# use mean of each sample??
# NO OB25 NORM DATA...
# choose O in PO instead
# only look at intersect genes
normdir<-"D:\\xnm\\work\\DATA\\NORM"
WJA<-"WJ_Aorta_NORM.txt"
HNA<-"HN_Aorta_NORM.txt"
PO<-"PO_NORM.txt"

WJA_dat<-read.table(paste0(normdir,"\\",WJA),header=T)
HNA_dat<-read.table(paste0(normdir,"\\",HNA),header=T)
PO_dat<-read.table(paste0(normdir,"\\",PO),header=T)

# 0.2.1 Without O, no need to aggrdup and merge
W<-normMean(WJA_dat,2,5,log2=T,aggr=F)
J<-normMean(WJA_dat,6,9,T,F)
H<-normMean(HNA_dat,2,5,T,F)
N<-normMean(HNA_dat,6,9,T,F)
WJHN<-cbind(W,J[,2],H[,2],N[,2])
colnames(WJHN)<-c("ID","W","J","H","N")
write.table(WJHN,"WJHN_log2_norm_mean_no_aggr.txt",quote = F,sep = "\t")


WJHN<-read.table("WJHN_log2_norm_mean_no_aggr.txt",sep = "\t")

corrPlot(WJHN,"W","J")
corrPlot(WJHN,"W","H")
corrPlot(WJHN,"H","N")

# to find out genes with high expression in the diagnal
test<-WJHN[order(WJHN$W,WJHN$H,decreasing = T),]
# From the plot, we find that expr>15 can be labeled


# 0.2.2 With O,  aggrdup and merge
aW<-normMean(WJA_dat,2,5,log2=T,aggr=T)
aJ<-normMean(WJA_dat,6,9,T,T)
aH<-normMean(HNA_dat,2,5,T,T)
aN<-normMean(HNA_dat,6,9,T,T)
P<-normMean(PO_dat,2,4,T,T)
O<-normMean(PO_dat,5,7,T,T)
aWJHN<-cbind(aW,aJ[,2],aH[,2],aN[,2]) # it is okay as the order of ID wont change after aggrDup()
colnames(aWJHN)<-c("ID","W","J","H","N")
write.table(aWJHN,"WJHN_log2_norm_mean_aggr.txt",quote = F,sep = "\t")
aWJHN<-read.table("WJHN_log2_norm_mean_aggr.txt",sep = "\t")

PO<-cbind(P,O[,2])
colnames(PO)<-c("ID","P","O")
dat_merge<-merge(aWJHN, PO, by="ID", all=F)
colnames(dat_merge)<-c("ID","W","J","H","N","P","O")
write.table(dat_merge,"WJHNO_log2_norm_mean_aggr.txt",quote = F,sep = "\t")

dat_merge<-read.table("WJHNO_log2_norm_mean_aggr.txt",sep = "\t")

# corrPlot(dat_merge,"W","H")
corrPlot(dat_merge,"W","O")
corrPlot(dat_merge,"H","O")
# corrPlot(dat_merge,"P","O")
corrPlot(PO,"P","O")

