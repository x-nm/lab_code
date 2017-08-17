# modified the sample code of clusterProfiler, for the server to use
# 2017.8.11 by xnm

## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

## ----echo=FALSE, results='hide', message=FALSE---------------------------
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase) # XNL problem in server
library(clusterProfiler)

## ------------------------------------------------------------------------
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
#eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db") # a map function for gene symbol
#'select()' returned 1:1 mapping between keys and columns
head(eg)

# SYMBOL ENTREZID
# 1   GPX3     2878
# 2   GLRX     2745
# 3    LBP     3929
# 4  CRYAB     1410
# 5  DEFB1     1672
# 6  HCLS1     3059


## ------------------------------------------------------------------------
library(org.Hs.eg.db) 
keytypes(org.Hs.eg.db) # seems like an ID database, can transform gene ID from different types. Seems useful.

# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"
# [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"
# [26] "UNIPROT"

## ------------------------------------------------------------------------
#ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), annoDb="org.Hs.eg.db")
head(ids)

# SYMBOL    UNIPROT         ENSEMBL
# 1   GPX3     P22352 ENSG00000211445
# 2   GLRX A0A024RAM2 ENSG00000173221
# 3   GLRX     P35754 ENSG00000173221
# 4    LBP     P18428 ENSG00000129988
# 5    LBP     Q8TCF0 ENSG00000129988
# 6  CRYAB     P02511 ENSG00000109846


## ------------------------------------------------------------------------
data(gcSample)
hg <- gcSample[[1]]
head(hg)

# [1] "4597"  "7111"  "5266"  "2175"  "755"   "23046"


eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
# bitr_kegg: not in the server version
head(eg2np)

## ----eval=FALSE----------------------------------------------------------
#  bitr_kegg("Z5100", fromType="kegg", toType='ncbi-geneid', organism='ece')

## ------------------------------------------------------------------------
bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

## ----warning=FALSE-------------------------------------------------------
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
# ENTREZID    UNIPROT         ENSEMBL
# 1     4312     B4DN15 ENSG00000196611
# 2     4312     P03956 ENSG00000196611
# 3     4312     Q53G95 ENSG00000196611
# 4     8318     O75419 ENSG00000093009
# 5    10874 A0A0B4J202 ENSG00000109255
# 6    10874     P48645 ENSG00000109255

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
# PROBLEM WITH OrgDb
head(ggo)

## ------------------------------------------------------------------------
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                # annoDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# GO TERMS enrichment, seems useful. But not really convenient to use.

head(ego)
# cannot head

# > ego
# 207 human Genes to  CC  test for over-representation.
# with p value < 0.01


ego <- enrichGO(gene          = x,
                universe      = names(geneList),
                # annoDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)


## ----eval=FALSE----------------------------------------------------------
#  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
#                  OrgDb         = org.Hs.eg.db,
#  		keytype       = 'ENSEMBL',
#                  ont           = "CC",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
#  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

## ----eval=FALSE----------------------------------------------------------
#  ego3 <- gseGO(geneList     = geneList,
#                OrgDb        = org.Hs.eg.db,
#                ont          = "CC",
#                nPerm        = 1000,
#                minGSSize    = 100,
#                maxGSSize    = 500,
#                pvalueCutoff = 0.05,
#                verbose      = FALSE)

## ------------------------------------------------------------------------
search_kegg_organism('ece', by='kegg_code') 
# no such function
ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)
head(ecoli)

## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

# > kk
# 207 hsa Genes to  KEGG  test for over-representation.
# with p value < 0.05

head(kk)
# not able to head()

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

# > kk2
# [1] "GSEA analysis result Object..."

head(kk2)
# not able to head()

## ----eval = FALSE--------------------------------------------------------
#  mkk <- enrichMKEGG(gene = gene,
#                     organism = 'hsa')

## ----eval=FALSE----------------------------------------------------------
#  mkk2 <- gseMKEGG(geneList = geneList,
#                   species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
#  david <- enrichDAVID(gene = gene,
#                       idType = "ENTREZ_GENE_ID",
#                       listType = "Gene",
#                       annotation = "KEGG_PATHWAY",
#                       david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile) # no such function

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
pdf("barplot_ggo.pdf")
barplot(ggo, drop=TRUE, showCategory=12)
dev.off()

## ----fig.height=5, fig.width=8-------------------------------------------
pdf("barplot_ego.pdf")
barplot(ego, showCategory=8)
dev.off()

## ------------------------------------------------------------------------
pdf("dotplot_ego.pdf")
dotplot(ego)
# Error in (function (classes, fdef, mtable)  :
#             unable to find an inherited method for function ¡®dotplot¡¯ for signature ¡®"logical"¡¯
#           
dev.off()

## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
#  enrichMap(ego)

## ----fig.height=14, fig.width=14, eval=FALSE-----------------------------
#  ## categorySize can be scaled by 'pvalue' or 'geneNum'
#  cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
pdf("gograph_ego.pdf")
plotGOgraph(ego)
dev.off()

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
pdf("gseaplot.pdf")
gseaplot(kk2, geneSetID = "hsa04145")
dev.off()

## ----eval=FALSE----------------------------------------------------------
#  browseKEGG(kk, 'hsa04110')

## ----eval=FALSE----------------------------------------------------------
#  library("pathview")
#  hsa04110 <- pathview(gene.data  = geneList,
#                       pathway.id = "hsa04110",
#                       species    = "hsa",
#                       limit      = list(gene=max(abs(geneList)), cpd=1))

## ------------------------------------------------------------------------
data(gcSample)
lapply(gcSample, head)

## ------------------------------------------------------------------------
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck)) # error

## ------------------------------------------------------------------------
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

# Entrez       FC       group othergroup
# 4312    4312 4.572613 upregulated          B
# 8318    8318 4.514594 upregulated          B
# 10874  10874 4.418218 upregulated          B
# 55143  55143 4.144075 upregulated          B
# 55388  55388 3.876258 upregulated          B
# 991      991 3.677857 upregulated          B


formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res)) # error

## ----fig.height=7, fig.width=9-------------------------------------------
pdf("dotplot_ck.pdf")
dotplot(ck)
dev.off()

## ----fig.height=6, fig.width=10------------------------------------------
pdf("dotplot_formula_res.pdf")
dotplot(formula_res)
dev.off()

pdf("dotplot_formula_res_plus.pdf")
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

# Error in .local(object, ...) : unused argument (x = ~group)

dev.off()

## ----echo=FALSE----------------------------------------------------------
sessionInfo()
