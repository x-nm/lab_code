# other figures
# 2017.07.07 by xnm

setwd("D:\\xnm\\work\\17.5.25\\LFC")
saveDir<-"D:\\xnm\\work\\17.5.25\\Small Figures/"
setwd(saveDir)

source("D:\\xnm\\work\\CODE\\figures_func.R")
options(stringsAsFactors = FALSE) # must, or ID will be transformed into num...

fdr<-0.1
fc<-0.5

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# 0.3 logFC的散点图（两两，区分>2什么的颜色，也可以不用）[ALL的话，就用原始数据吧]
# DEG a/b/both情况进行着色。阈值，形状
# （展现观察组间的相似性或差异性）【may also draw a separate one for DEGs】
# 给找出的基因标颜色。另外，右上方区域内的基因，不妨看一看
W_dat<-read.table("WJA_LFC.txt",sep = "\t")
H_dat<-read.table("HNA_LFC.txt",sep = "\t")
OB25_dat<-read.table("OB25_LFC.txt",sep = "\t")
PO_dat<-read.table("PO_LFC.txt",sep = "\t")


W_dat<-aggrDup(W_dat)
H_dat<-aggrDup(H_dat)
OB25_dat<-aggrDup(OB25_dat)

WH_dat<-merge(W_dat,H_dat,by="ID",all=F)
WH_dat<-WH_dat[apply(WH_dat[,-1],1,function(x) !all(is.na(x))),]
colnames(WH_dat)<-c("ID","WJ","WJ.p-val","HN","HN.p-val")

dat<-merge(WH_dat,OB25_dat,by="ID",all=F)
dat<-merge(dat,PO_dat,by="ID",all=F)
colnames(dat)<-c("ID","WJ","WJ.p-val","HN","HN.p-val","OB25","OB25.p-val","PO","PO.p-val")

setwd(saveDir)
write.table(WH_dat, "WH_lfc_aggr.txt", quote = F, sep = "\t")
write.table(dat, "WHPOB25_lfc_aggr.txt", quote = F, sep = "\t")

####

WH_dat<-read.table("WH_lfc_aggr.txt",sep = "\t")
dat<-read.table("WHPOB25_lfc_aggr.txt",sep = "\t")



# where is r-square?
# na --> no r-square

corrPlot(WH_dat,"WJ","HN","lfc")
corrPlot(dat,"WJ","OB25","lfc")
corrPlot(dat,"HN","OB25","lfc")
corrPlot(dat,"WJ","PO","lfc")
corrPlot(dat,"HN","PO","lfc")
corrPlot(dat,"OB25","PO","lfc")



# DEG_LFC
W_DEG<-dat_DEG(W_dat,fdr = fdr,fc = fc)
H_DEG<-dat_DEG(H_dat,fdr = fdr,fc = fc)
OB25_DEG<-dat_DEG(OB25_dat,fdr = fdr,fc = fc)
PO_DEG<-dat_DEG(PO_dat,fdr,fc)

WH_deg<-merge(W_DEG,H_DEG,by="ID",all=TRUE)
WH_deg<-WH_deg[apply(WH_deg[,-1],1,function(x) !all(is.na(x))),]
colnames(WH_deg)<-c("ID","WJ","WJ.p-val","HN","HN.p-val")

dat_deg<-merge(WH_deg,OB25_DEG,by="ID",all=TRUE)
dat_deg<-merge(dat_deg,PO_DEG,by="ID",all=TRUE)
colnames(dat_deg)<-c("ID","WJ","WJ.p-val","HN","HN.p-val","OB25","OB25.p-val","PO","PO.p-val")

setwd(saveDir)
write.table(WH_deg, "WH_lfc_deg_aggr.txt", quote = F, sep = "\t")
write.table(dat_deg, "WHPOB25_lfc_deg_aggr.txt", quote = F, sep = "\t")


WH_deg<-read.table("WH_lfc_deg_aggr.txt",sep = "\t")

corrPlot(WH_deg,"WJ","HN","lfc_deg")
corrPlot(dat_deg,"WJ","OB25","lfc_deg")
corrPlot(dat_deg,"HN","OB25","lfc_deg")
corrPlot(dat_deg,"HN","PO","lfc_deg")
corrPlot(dat_deg,"WJ","PO","lfc_deg")
corrPlot(dat_deg,"OB25","PO","lfc_deg")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# 0.4 每个数据组的logfc和-log10pval（y）的图
# （在0和1.3的位置画一条虚线【看重要的基因】【may also draw a separate one for DEGs】
# WH_dat, WH_deg, dat, dat_deg: W,H,P,OB25


lfc_lgp(WH_dat,"WJ","all")
lfc_lgp(WH_deg,"WJ","DEG")
lfc_lgp(WH_dat,"HN","all")
lfc_lgp(WH_deg,"HN","DEG")
lfc_lgp(dat,"PO","all")
lfc_lgp(dat_deg,"PO","DEG")
lfc_lgp(dat,"OB25","all")
lfc_lgp(dat_deg,"OB25","DEG")







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# 0.5 VENN of DEG intersections of W-H-OB25
dat_all<-wholeDataSet(fdr,fc,samples, aggr = T)
W<-subset(dat_all,subset = dat_all$sample=="WJA")$ID
H<-subset(dat_all,subset = dat_all$sample=="HNA")$ID
OB25<-subset(dat_all,subset = dat_all$sample=="OB25")$ID
W<-na.omit(W)
H<-na.omit(H)

T<-venn.diagram(list(W=W,H=H,OB25=OB25),
                filename=NULL ,imagetype = "png",
                # print.mode = "percent", # or default "raw"
                # main="DEG Venn of W-H-OB25",
                lwd=1,lty=2,
                col=c('#abd9e9','#fdae61','#2c7bb6'),fill=c('#abd9e9','#fdae61','#2c7bb6') ,
                category.names = c("","",""),
                rotation.degree=180) 
grid.draw(T)

# paste0(saveDir,"venn-WHOB25.png")


# hit gene list
a<-W[(W%in%H)]
b<-W[W%in%OB25]
c<-H[H%in%OB25]

hit3<-a[a%in%b]
hit2<-c(a,b,c)
hit23<-hit2[!duplicated(hit2)]
hit2<-hit2[!hit2%in%hit3]

a<-cbind(a,rep("WH",length(a)))
b<-cbind(b,rep("WOB",length(b)))
c<-cbind(c,rep("HOB",length(c)))
hit2_lb<-as.data.frame(rbind(a,b,c))
colnames(hit2_lb)<-c("ID","HITS_LABLE")
hit2_lb<-hit2_lb[!duplicated(hit2_lb$ID),]
hit2_lb<-hit2_lb[!hit2_lb$ID%in%hit3,]
save(hit2,hit23,hit3,hit2_lb,file = paste0(saveDir,"hit23.Rdata"))


load(paste0(saveDir,"hit23.Rdata"))

write(hit3,paste0(saveDir,"WHOB25-hit3.txt"))
write(hit23,paste0(saveDir,"WHOB25-hit23.txt"))
write(hit2,paste0(saveDir,"WHOB25-hit2.txt"))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 0.6 hit3 colored LFC-corr of WH

dat_new<-WH_dat[,c("ID","WJ","HN")]
index<-dat_new$ID%in%hit3
dat_new<-cbind(dat_new,index)
colnames(dat_new)[4]<-"DEG_of_WHOB25"



load("dat_new_for_fig6.Rdata")

title <- paste0("Correlation of LogFC of WH")
ggplot(dat_new, aes(x = WJ, y = HN, color=DEG_of_WHOB25))+
  geom_point()+ # or no shape = 1
  geom_smooth(method=lm)+ # or indianred3
  labs(x="W_LFC", y="H_LFC" )+
  ggtitle(title)+
  theme(
    plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
    axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
  )
# annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))
ggsave(paste0("WH_lfc_corr_colored.png"),width = 5.5, height = 4)

dat_deg<-WH_deg[,c("ID","WJ","HN")]
index<-dat_deg$ID%in%hit3
dat_deg<-cbind(dat_deg,index)
colnames(dat_deg)[4]<-"DEG_of_WHOB25"
save(dat_new,dat_deg,file="dat_new_for_fig6.Rdata")




title <- paste0("Correlation of LogFC of DEGs of WH")
ggplot(dat_deg, aes(x = WJ, y = HN, color=DEG_of_WHOB25))+
  geom_point()+ # or no shape = 1
  geom_smooth(method=lm)+ # or indianred3
  labs(x="W_LFC", y="H_LFC" )+
  ggtitle(title)+
  theme(
    plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
    axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
  )
  # annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))
ggsave(paste0("WH_lfc_deg_corr_colored.png"),width = 5.5, height = 4)
