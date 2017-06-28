# t-test for edge fts
# 2016.11.6

# save(matExpr, matCtl,geneID, file=paste("edgeFts_",label,".RData", sep = ""))

################################################
# setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\WCGNA\\Data")
setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DE")

label <- "HN"
FILE <- paste("edgeFts_",label,".RData", sep = "")

lname <- load(FILE)
lname
################################################

ftNum<-dim(matExpr)[1]

pVal<-c()
for(i in 1:ftNum){
  t<-t.test(matExpr[i,],matCtl[i,])
  pVal<-rbind(pVal,t$p.value)
}

index<-pVal<0.05

# matrix after t-test
matTTestExpr<-matExpr[index,]
matTTestCtl<-matCtl[index,]

# select and annotate
uvSel<-matExpr[index,1:2]
uvSelNum<-dim(uvSel)[1]

u2geneID=geneID[uvSel[,1]]
v2geneID=geneID[uvSel[,2]]

annot = read.csv("gene_id2symbol.txt", header = F, sep = "\t")
names(annot) <- c("geneID","geneSymbol")
u2geneSym = match(u2geneID, annot$geneID)
v2geneSym = match(v2geneID, annot$geneID)

uvSelSym <- data.frame(edge=row.names(uvSel),
                       uGeneSym = annot$geneSymbol[u2geneSym],
                       vGeneSym = annot$geneSymbol[v2geneSym])

save(matTTestExpr,matTTestCtl,uvSelSym,file=paste("DETTest_",label,".RData",sep=""))

