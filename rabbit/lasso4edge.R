# lasso for edge
# TRAIL VERSION
# 2016.11.6

# save(matTTestExpr,matTTestCtl,uvSelSym,file=paste("DETTest_",label,".RData",sep=""))

################################################
# setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\WCGNA\\Data")
setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DE")
library(glmnet)

label <- "HN"
FILE <- paste("DETTest_",label,".RData",sep="")

lname <- load(FILE)
lname
################################################

# outcome <- read.csv("LASSOoutcome.csv") # a col of pheno
pheno <- as.data.frame(c(rep(1,4),rep(0,4)))
colnames(pheno)<-"pheno"
# gene <- read.csv("LASSOcovariate.csv") # expr file, sp#row, gene#col
# head(outcome)
gene <- t(cbind(matTTestExpr[,3:6],matTTestCtl[,3:6]))
# is.matrix(gene)
# dim(gene)  # n=99, P=28
X <- as.matrix(gene)
# attach(pheno)
# outcome<-pheno
# attach(outcome)
pheno<-as.matrix(pheno)

set.seed(1234)
# model <- cv.glmnet(X, group, family="binomial", type.measure="deviance")
model <- cv.glmnet(X, pheno, family="binomial", type.measure="deviance")
plot(model)  # CV plot for lambda
model$lambda.1se # find optimal lambda

model.final <- model$glmnet.fit
model.coef <- coef(model$glmnet.fit, s = model$lambda.1se)
all.coef <- coef(model$glmnet.fit, s =  min(model.final$lambda))
index <- which(model.coef != 0)
cbind(rownames(model.coef)[index], round(model.coef[index],4)) #选出的features


edgeLasso<-uvSelSym[index-1,]

save(edgeLasso,file = paste("edgeLasso_",label,".RData",sep = ""))
