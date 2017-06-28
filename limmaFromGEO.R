# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Mar 7 06:08:44 EST 2017
# MODIFIED FOR LOCAL USE BY XNM IN 2017.3.8

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#                                DEFINATION                                   #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
setwd("D:\\xnm\\work\\DATA")
gset <- getGEO("GSE37558", GSEMatrix =TRUE, AnnotGPL=TRUE)
DATA<-gset


# label<-"OB12"
# gsms <- "0000111XXXXXXXXXXXXXXXXXXXXXXXXX" #OB2
# gsms <- "0000XXX111XXXXXXXXXXXXXXXXXXXXXX" #OB8
# gsms <- "0000XXXXXX111XXXXXXXXXXXXXXXXXXX" #OB12
# gsms <- "0000XXXXXXXXX111XXXXXXXXXXXXXXXX" #OB25

label<-"CV25"
# gsms <- "XXXXXXXXXXXXXXXX0000111XXXXXXXXX" #CV2
# gsms <- "XXXXXXXXXXXXXXXX0000XXX111XXXXXX" #CV8
# gsms <- "XXXXXXXXXXXXXXXX0000XXXXXX111XXX" #CV12
gsms <- "XXXXXXXXXXXXXXXX0000XXXXXXXXX111" #CV25




# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("GEOquery")

# load series and platform data from GEO

# gset <- getGEO("GSE37558", GSEMatrix =TRUE, AnnotGPL=TRUE)
gset<-DATA
if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1 # GPL: platform
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
# gsms <- "0000111XXXXXXXXXXXXXXXXXXXXXXXXX" #OB2
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# save results
num<-nrow(ex)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=num)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","GO.Function","GO.Function.ID","GO.Process","GO.Process.ID" ))
write.table(tT, paste(label,"lfc_raw", sep = ""), row.names=F, sep="\t", quote = F)


# save LFC data, also do some cleaning
lfc<-subset(tT,select=c("Gene.symbol","logFC","adj.P.Val"))
colnames(lfc)<-c("ID","LFC","padj")

# count item numbers under different padj thresholds
padjCt<-function(data){
  return(c(nrow(subset(data, data$pdaj<0.01)),
           nrow(subset(data, data$pdaj<0.05)),
           nrow(subset(data, data$pdaj<0.1))))
}

count<-c()
count<-rbind(count,c("state","<0.01","<0.05","<0.1"))
count<-rbind(count,c("Before",padjCt(lfc)))

# clean: del NA ID 
emt<-which(lfc$ID=="")
lfc<-lfc[-emt,]

# clean: deal with duplicated ID
dup<-which(duplicated(lfc$ID))
dupItems<-lfc[dup,]
# lfc_m: just minus the duplicated items 
lfc_m<-lfc[-dup,]
# lfc_a: add up duplicted items' value
lfc_a<-lfc_m
mch<-match(dupItems$ID,lfc_m$ID)
for(i in 1:length(dup)){
  a<-lfc_m[mch[i],]$LFC
  b<-dupItems[i,]$LFC
  if(!is.na(a)){
    if(!is.na(b)) lfc_a[mch[i],]$LFC<-a+b
  }
  else{
    if(!is.na(b)) lfc_a[mch[i],]$LFC<-b
  }
}
# SHOULD CARE ABOUT PADJ?
count<-rbind(count,c("After",padjCt(lfc_m)))

write.table(lfc_m,paste(label,"m_LFC.txt", sep = ""),quote = F, row.names = F, sep = "\t")
write.table(lfc_a,paste(label,"a_LFC.txt", sep = ""),quote = F, row.names = F, sep = "\t")
write.table(count,paste(label,"count.txt", sep = ""),quote = F, row.names = F, sep = "\t", col.names = F)
