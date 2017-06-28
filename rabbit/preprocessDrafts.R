
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)


# ��deseq2��ķ���ȥ��count<=1�Ļ���
nrow(dds)
# [1] 41559
ddsShrunk <- dds[rowSums(counts(dds))>1,]
dim(ddsShrunk)
# [1] 28945     8


#####################
dds <- DESeq(dds)
res <- results(dds)

vsd<-varianceStabilizingTransformation(dds)
nm<-assay(vsd)
head(nm)
#head(counts(dds))
rld<-rlog(dds)
nmRld<-assay(rld)
head(nmRld)

####################
ddsShr<-DESeq(ddsShrunk)
vsdShr<-varianceStabilizingTransformation(ddsShr)
nmVsrShr<-assay(vsdShr)
head(nmVsrShr)
head(counts(ddsShr))

####################
ddsShrunk2 <- dds[rowSums(counts(dds))>8,]
dim(ddsShrunk2)
# [1] 23819     8

####################
head(counts(dds, normalized = T))
head(counts(dds, normalized = F))


####################
cds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds, normalized = T))
head(counts(cds, normalized = F))