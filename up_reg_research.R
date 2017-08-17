# upstream regulator research
# 2017.7.20 by xnm

setwd("D:\\xnm\\work\\DATA\\IPA")

# table of top25
WJ<-read.table("WJ_upstream.txt",sep="\t")
WJ_top25<-WJ[order(WJ$Activation.z.score,decreasing = T),][1:25,1:6]
write.table(WJ_top25,"WJ_top25.txt",sep = "\t",quote = F,row.names = F)
HN<-read.table("HN_upstream.txt",sep="\t")
HN_top25<-HN[order(HN$Activation.z.score,decreasing = T),][1:25,1:6]
write.table(HN_top25,"HN_top25.txt",sep = "\t",quote = F,row.names = F)

# heatmap of top25
hm_dat1<-cbind(WJ_top25[,c(1,5)],rep("WJ",25))
colnames(hm_dat1)[3]<-"Sample"
hm_dat2<-cbind(HN_top25[,c(1,5)],rep("HN",25))
colnames(hm_dat2)[3]<-"Sample"
hm_dat<-rbind(hm_dat1,hm_dat2)  

ID<-rev(na.omit(hm_dat$Upstream.Regulator))

hm_dat<-transform(hm_dat,Upstream.Regulator = factor(Upstream.Regulator,levels = unique(ID)))

title<-"Top 25 Upstream Regulators"

library(ggplot2)
base_size<-11
hColor<-"red"
mColor<-"white"
lColor<-"steelblue2"

ggplot(hm_dat,aes(x=Sample,y=Upstream.Regulator))+
  geom_tile(aes(fill=Activation.z.score))+
  scale_fill_gradient2(mid=mColor, high=hColor, low=lColor) +
  xlab("")+
  ylab("")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 1.2*base_size, hjust = 0.5,face = "bold"),
    axis.text.x = element_text(angle = 90,  vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
    axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50"),
    legend.title = element_text(size = 0.8 *  base_size),
    legend.text = element_text(size = 0.8 *  base_size),
    legend.key.size = unit(0.8, "cm")
    # axis.line.x = element_line(color = level$starLevel,linetype = 2, size = 0.2 * base_size)
  )+
  ggtitle(title)+ # title
  coord_fixed()
