
library(VennDiagram)
setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DEG\\heatmap\\input")

KG <- read.table("keyGeneList.txt")$V1
BM <- read.table("BM.txt")$V1
posBM <- read.table("posBM.txt")$V1
negBM <- read.table("negBM.txt")$V1
OP <- read.table("osporoList.txt")$V1
HM <- read.table("hyclMutList.txt")$V1
WM <- read.table("WHHLMutList.txt")$V1
obd <- read.table("ostbDiff.txt")$V1
ctd <- read.table("ctlgDv.txt")$V1
osf <- read.table("osfct.txt")$V1
LP <- read.table("lipidMetList.txt")$V1
IN <-read.table("inflmResList.txt")$V1

setwd("D:\\Files\\Study\\交大\\课题\\16.10.27\\DEG\\heatmap\\output")

LIST <- list(Keygenes=KG, 
             #BM=BM, 
             #posBM=posBM, 
             #negBM=negBM#, 
             #Osporosis=OP, 
             #HM=HM, 
             #WM=WM#, 
             # OsteoblastDiff=obd, 
             # CartilageDev = ctd,
             # Ossification = osf,
             # LP=LP,
             IN = IN
             )




T<-venn.diagram(LIST, filename = NULL, lwd=1,lty=2 )#,
                #col=c('red','green')
                #fill=c('red','green')
                #cat.col=c('red','green') ,
                #rotation.degree=90) 
grid.draw(T) 

