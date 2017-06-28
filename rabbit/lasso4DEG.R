# lasso for DEG
# TRAIL VERSION
# 2016.11.6

#############################################################################
# save(NORM1, NORM2, DEG1, DEG2, file="normValOfDEG.RData")

setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DEG\\lasso")
lname<-load("normValOfDEG.RData")

library(glmnet)

#############################################################################

# #############################################################################
# # GET DEG EXPR DATA
# setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DEG\\heatmap\\input")
# 
# DEG1 <- read.csv("DEG_WJ.txt",header = T,sep = "\t")# WJ as 1, HN as 2
# DEG1 <- data.frame(geneID=row.names(DEG1),geneSymbol=DEG1[,7]) # GENE ID AS ROW NAME, GENE SYMBOL AS CONTENT
# DEG2 <- read.csv("DEG_HN.txt",header = T,sep = "\t")
# DEG2 <- data.frame(geneID=row.names(DEG2),geneSymbol=DEG2[,7])
# 
# NORM1 <- read.csv("normValue_WJ.txt", header = T, sep = "\t")
# NORM2 <- read.csv("normValue_HN.txt",header = T, sep = "\t")
# 
# indexInNorm1 <- match(DEG1$geneID,row.names(NORM1))
# indexInNorm2 <- match(DEG2$geneID,row.names(NORM2))
# 
# NORM1<-NORM1[indexInNorm1,]
# sum(is.na(NORM1)) #2264/8=283, 283 DEGs with low counts
# NORM1[is.na(NORM1)]<-0
# NORM2<-NORM2[indexInNorm2,]
# sum(is.na(NORM2)) #3624/8=453
# NORM2[is.na(NORM2)]<-0
# 
# setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DEG\\lasso")
# save(NORM1, NORM2, DEG1, DEG2, file="normValOfDEG.RData")
# 
# ##########################################################################

lasso<-function(X,pheno){
  set.seed(1234)
  model <- cv.glmnet(X, pheno, family="binomial", type.measure="deviance")
  plot(model)  # CV plot for lambda
  model$lambda.1se # find optimal lambda
  
  model.final <- model$glmnet.fit
  model.coef <- coef(model$glmnet.fit, s = model$lambda.1se)
  all.coef <- coef(model$glmnet.fit, s =  min(model.final$lambda))
  index <- which(model.coef != 0)
  cbind(rownames(model.coef)[index], round(model.coef[index],4)) #选出的features
  return(index)
}

# LASSO
pheno <- as.data.frame(c(rep(1,4),rep(0,4)))
colnames(pheno)<-"pheno"
pheno<-as.matrix(pheno)

X1 <- as.matrix(t(NORM1))
X2 <- as.matrix(t(NORM2))

index1<-lasso(X1,pheno)
index2<-lasso(X2,pheno)

DEGLasso1<-DEG1[index1-1,]
DEGLasso2<-DEG2[index2-1,]

save(DEGLasso1, DEGLasso2,file = "DEGLasso.RData")
