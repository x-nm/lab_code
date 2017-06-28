# 2016.10.27 BY xnm
# INPUT: "calPCC_HN.RData"

# 1. CAL PCC, IN calPCC_linux.R
# save(geneNum, spNumCtl, spNumExpr, dataExpr, dataCtl, uvPairs, file = paste("calPCC_",label, ".RData",sep=""))


################################################
# setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\WCGNA\\Data")
setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DE")

label <- "HN"
FILE <- paste("calPCC_",label, ".RData",sep="")

lname <- load(FILE)
lname
################################################

# 2. Construct edge feature matrix

# DATA:geneNum ROW, spNum COL. 

uvNum <- dim(uvPairs)[1]
uvPairs<-as.data.frame(uvPairs)

matExpr <- c()
matCtl <- c()

K1 <- sqrt((spNumCtl-1)/spNumCtl)
K2 <- sqrt((spNumExpr-1)/spNumExpr)

for (i in 1:uvNum){
  u<- uvPairs[i,1]
  v<- uvPairs[i,2]
  uExpr <- as.numeric(dataExpr[u,])
  uCtl <- as.numeric(dataCtl[u,])
  vExpr <- as.numeric(dataExpr[v,])
  vCtl <- as.numeric(dataCtl[v,])
  uEM <- mean(uExpr, na.rm = T)
  uCM <- mean(uCtl, na.rm = T)
  vEM <- mean(vExpr, na.rm = T)
  vCM <- mean(vCtl, na.rm = T)
  uESd <- sd(uExpr, na.rm = T)
  uCSd <- sd(uCtl, na.rm = T)
  vESd <- sd(vExpr, na.rm = T)
  vCSd <- sd(vCtl, na.rm = T)
  
  ftE4E <- ((uExpr - uEM)*(vExpr - vEM))/(K1 * K1 * uESd * vESd) #feature value for Expr, cal from expr formula
  ftC4E <- ((uExpr - uCM)*(vExpr - vCM))/(K2 * K2 * uCSd * vCSd) #feature value for Expr, cal from Ctl formula
  ftE4E <- rbind(u,v, as.matrix(ftE4E)) # use rbind rather than cbind, as vector will be converted as a column here.
  ftC4E <- rbind(u,v, as.matrix(ftC4E))
  colnames(ftE4E)<-paste("E",u,"-",v,sep = "")
  colnames(ftC4E)<-paste("C",u,"-",v,sep = "")
  matExpr <- cbind(matExpr, ftE4E, ftC4E)
  
  ftE4C <- ((uCtl - uEM)*(vCtl - vEM))/(K1 * K1 * uESd * vESd)
  ftC4C<- ((uCtl - uCM)*(vCtl - vCM))/(K2 * K2 * uCSd * vCSd)
  ftE4C <- rbind(u,v, as.matrix(ftE4C))
  ftC4C <- rbind(u,v, as.matrix(ftC4C))
  colnames(ftE4C)<-paste("E",u,"-",v,sep = "")
  colnames(ftC4C)<-paste("C",u,"-",v,sep = "")
  matCtl <- cbind(matCtl, ftE4C, ftC4C)
}

matCtl<-t(matCtl)
colnames(matCtl)<-c("u","v",colnames(dataCtl))
matExpr<-t(matExpr)
colnames(matExpr)<-c("u","v",colnames(dataExpr))

geneID<-row.names(dataExpr)

save(matExpr, matCtl,geneID, file=paste("edgeFts_",label,".RData", sep = ""))

