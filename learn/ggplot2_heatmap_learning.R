# 2017.05.25 by xnm


ggplot(m.WHO, aes(variable, Country)) +
  #��tile�����л�����ͼ
  geom_tile(aes(fill=rescale)) +
  scale_fill_gradient2(mid="white", high="red", low="blue", name = "LFC") +
  labs(x="Country", y="Index", face = "bold") +
  theme_bw() +
  theme(
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=12, colour="grey50"),
    axis.text.y=element_text(size=12, colour="grey50"),
    legend.title=element_text(size=14),
    legend.text=element_text(size=12),
    legend.key.size = unit(0.8, "cm"))#��Ҫ����grid��������legend�Ĵ�С


#########################################################################
setwd("D:\\xnm\\work\\17.5.25\\LFC/")

dat1<-read.csv(file = "HN_Aorta_LFC.txt", header = T, sep = "\t")
dat2<-read.csv(file = "WJ_Aorta_LFC.txt", header = T, sep = "\t")
dat3<-read.csv("OB25_LFC.txt",sep = "\t")
# dat1<-arrange(dat1, desc(ID))
# dat2<-arrange(dat2, desc(ID))
# dat3<-arrange(dat3, desc(ID))
dat<-rbind(dat1[1:300,],dat2[1:300,],dat3[1:300,])
dat<-cbind(dat,c(rep("HNA",300),rep("WJA",300),rep("HNH",300)))
colnames(dat)[4]<-"sample"


dat<-transform(dat, sample=factor(sample, levels = unique(sample)))
dat<-transform(dat, ID=factor(ID, levels = unique(ID)))


plot<-ggplot(dat, aes(sample, ID)) +
  geom_tile(aes(fill=LFC)) +
  scale_fill_gradient2(mid="white", high="red", low="green", name = "LFC") +
  labs(x="samples", y="Gene", face = "bold") +
  theme_bw() +
  theme(
    axis.title.x=element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12, color = "grey50"),
    axis.text.y = element_text(size = 10, colour = "grey50"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"))

plot
ggsave("plot.pdf", plot, width = 3.15, height = 3.15) 


########### Learn from sample 2 ####################

# ���ಢ���򣬲�������������֯����
hc<-hclust(dist(data))
rowInd<-hc$order
hc<-hclust(dist(t(data)))
colInd<-hc$order
data.m<-data[rowInd,colInd]

base_size<-12

# main part, the same with former sample
p <- ggplot(dat, aes(sample, ID)) + 
  geom_tile(aes(fill=LFC)) +
  scale_fill_gradient2(mid="white", high="red", low="steelblue2", name = "LFC")


p + 
  theme_grey(base_size = base_size) + 
  labs(x = "",y = "") + 
  #����xlabel��ylabelΪ��
  scale_x_continuous(expand = c(0, 0),labels=coln,breaks=1:length(coln)) + 
  #����x������չ����Ϊ0���̶�Ϊ֮ǰ����Ʒ��
  scale_y_continuous(expand = c(0, 0),labels=rown,breaks=1:length(rown)) + 
  #����y������չ����Ϊ0���̶�Ϊ֮ǰ�Ļ�����
  opts(axis.ticks = theme_blank(), 
       axis.text.x = theme_text(size = base_size * 0.8, angle = 90, hjust = 0, colour = "grey50"), 
       axis.text.y = theme_text(size = base_size * 0.8, hjust=1, colour="grey50"))
#������������Ϊ��׼��0.8������������Խڣ�x������ת90�ȣ�ɫ��Ϊ�л�