# for ora results
setwd("D:\\xnm\\work\\17.2.15\\KeyGene\\WHP")

hn<-read.table("HN_ORA.txt",header = T, sep = "\t")
wj<-read.table("WJ_ORA.txt",header = T, sep = "\t")
po<-read.table("PO_ORA.txt",header = T, sep = "\t")
colnames(hn)

sum(wj[,9]<0.05)
sum(hn[,9]<0.05)
sum(po[,9]<0.05)

sum(wj[,9]<0.01)
sum(hn[,9]<0.01)
sum(po[,9]<0.01)

sum(wj[,9]<0.1)
sum(hn[,9]<0.1)
sum(po[,9]<0.1)

bind<-cbind(wj[,c(1,4,5,8,9)],hn[,c(1,4,5,8,9)],po[,c(1,4,5,8,9)])


bind2<-cbind(wj[,1],wj[,4]/wj[,5],wj[,c(8,9)],hn[,1],hn[,4]/hn[,5],hn[,c(8,9)],po[,1],po[,4]/po[,5],po[,c(8,9)])
boxplot(wj[,4]/wj[,5])
boxplot(hn[,4]/hn[,5])
boxplot(po[,4]/po[,5])


bind3<-merge(wj[,c(1,4,9)],hn[,c(1,4,9)],by="GeneSet")
bind3<-merge(bind3,po[,c(1,4,5,9)],by="GeneSet")

write.table(bind3,"ora.bind.txt", quote = F, sep = "\t")
