# script2 of otherTissueDataProcess
# 2017.05.11 by xnm




#######################################################################################
# 2. DEG; GET sorted results AND norm value.
#######################################################################################

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

summ<-file(paste0("output/", query, "_summaries.txt"))
# rm(summ)
sink(summ, append = T)
print(paste0(">",query))
summary(res)
sink()
# write.table(cat(readLines("summaries.log"), sep="\n"),file=paste0(query, "_summaries.txt"))
rm(summ)

#######################################################################################


#######################################################################################
# 3. WJ_KIDNEY_LFC(WITH ALL COLUMNS), ID + GENE SYMBOL, input data for heatmap drawing
# ID  GENE  baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
#######################################################################################

geneId2Sym <- read.csv("gene_id2symbol.txt", header = F, sep = "\t")

resOrdered <- res[order(res$padj),]
resOrdered <- as.data.frame(resOrdered)
index<-match(row.names(resOrdered),geneId2Sym[,1]) #是row.names吗？似乎应该是哦=-=
ID<-geneId2Sym[index,2]
LFC<-cbind(ID[1:41554],resOrdered[1:41554,]) 

indexExclude<-which(LFC[,1]=="")
LFC<-LFC[-indexExclude,] # "" excluded, 22682 left (not all like that, left with different # of NAs)
colnames(LFC)[1]<-"ID"

write.table(LFC,file=paste("output/",query, "_LFC.txt", sep = ""), sep = "\t", quote = F)

#######################################################################################


#######################################################################################
# 4. WJ_AORTA_NORM, ID + GENE SYMBOL, input data for heatmap drawing
# ID  GENE  R1  R2  R3  R4
#######################################################################################

NORM <- counts(dds, normalized = T)
NORM <- cbind(geneId2Sym[,2],NORM[1:41554,])
indexExclude<-which(NORM[,1]=="")
NORM <- NORM[-indexExclude,] # NA excluded, 22682 left
colnames(NORM)[1] <- "ID"

write.table(NORM,file=paste("output/",query, "_NORM.txt", sep = ""), quote = F, sep = "\t")

#######################################################################################