# view the state of tissue cell info
# 2017.06.26 by xnm

setwd("D:\\xnm\\work\\17.5.25\\LifeMap")

bone<-read.csv("bone_cell_info.txt",sep = "\t")
cart<-read.csv("cartilage_cell_info.txt",sep = "\t")

sum(bone$isPrimaryProgenitor)
sum(bone$isMainprimaryProgenitor)
sum(cart$isPrimaryProgenitor)
sum(cart$isMainprimaryProgenitor)

boxplot(bone[,c(5:8)])
summary(bone[,c(5:8)])

boxplot(cart[,c(5:8)])
summary(cart[,c(5:8)])

s1<-bone$Name
s2<-cart$Name

i1<-s1%in%s2
sum(i1)
i2<-s2%in%s1
sum(i2)


inter<-bone[i1,]
inter$Tissue<-rep("both",nrow(inter))
b<-bone[!i1,]
b$Tissue<-rep("bone",nrow(b))
c<-cart[!i2,]
c$Tissue<-rep("cartilage",nrow(c))

whole<-rbind(inter,b,c)
write.table(whole,"boneCart_CellInfo.txt",sep = "\t", quote = F,row.names = F)
