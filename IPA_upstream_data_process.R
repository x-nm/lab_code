# 2017.4.19 by xnm

setwd("D:\\xnm\\work\\17.1.4\\IPA\\17.4.18\\WH")
options(stringsAsFactors = F)

wh<-read.table("WH-upstream_analysis.txt",header = T, sep = "\t")
wh$Analysis<-substr(wh$Analysis,1,2)

moleType<-levels(as.factor(wh$Molecule.Type))
# 27 moleTypes, select from 11 to 18, and 21~27, all proteins with a single gene
# complex(10) is special, should look at it separately

pr<-subset(wh,wh$Molecule.Type%in%moleType[c(11:18,21:27)])
# <-merge(wh,moleType[c(11:18,21:27)],by.x="Molecule.Type", by.y = 1)
WJ<-subset(pr,pr$Analysis=="WJ")[,2:8]
HN<-subset(pr,pr$Analysis=="HN")[,2:8]


# A<-strsplit(WJ$Target.molecules.in.dataset[1],split = ",")
# A<-as.vector(A[[1]])
# length(A)

write.table(WJ, "WJ_upstream.txt", quote = F, sep = "\t")
write.table(HN, "HN_upstream.txt", quote = F, sep = "\t")
