library(ggplot2)
hc<-hclust(dist(data))
rowInd<-hc$order
hc<-hclust(dist(t(data)))
colInd<-hc$order
data.m<-data[rowInd,colInd] #���������������Ϊ��ɫ�鼯�У���ʾЧ���á���������Ͷ���Ʒ�з��飬
# ���������򣬾Ϳ���������һ����

data.m<-apply(data.m,1,rescale) #����Ϊ��׼�����ݽ��б任��ʹÿһ�ж���ɣ�0,1��֮������֡�
# �任�ķ���������scale,rescale�ȵȣ������Լ�����Ҫ���任��

data.m<-t(data.m) #�任�Ժ�ת���ˡ�
coln<-colnames(data.m) 
rown<-rownames(data.m) #������Ʒ���������ơ���Ϊgeom_tile������ǰ��������ţ�
# ������Ҫʹ�����ְ����ǵ����й̶�������

colnames(data.m)<-1:ncol(data.m)
rownames(data.m)<-1:nrow(data.m)
data.m<-melt(data.m) #ת�����ݳ��ʺ�geom_tileʹ�õ���ʽ
head(data.m)
# X1 X2     value
# 1  1  1 0.1898007
# 2  2  1 0.6627467
# 3  3  1 0.5417057
# 4  4  1 0.4877054
# 5  5  1 0.5096474
# 6  6  1 0.2626248

base_size<-12 #����Ĭ�������С��������Ʒ���߻���Ķ��ٶ�΢�䡣


p <- ggplot(data.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  #�趨������Ϊ��ǰ���У�������Ϊ��ǰ���У����ɫΪת���������
  scale_fill_gradient(low = "white",high = "steelblue")
  #�趨����ɫ�ĵ�ֵΪ��ɫ����ֵΪ����ɫ��


p + 
  theme_grey(base_size = base_size) + 
  labs(x = "", #����xlabel��ylabelΪ��
       y = "") + 
  scale_x_continuous(expand = c(0, 0),labels=coln,breaks=1:length(coln)) + 
  #����x������չ����Ϊ0���̶�Ϊ֮ǰ����Ʒ��
  scale_y_continuous(expand = c(0, 0),labels=rown,breaks=1:length(rown)) + 
  #����y������չ����Ϊ0���̶�Ϊ֮ǰ�Ļ�����
  opts(axis.ticks = theme_blank(), 
       axis.text.x = theme_text(size = base_size * 0.8, angle = 90, hjust = 0, colour = "grey50"), 
       axis.text.y = theme_text(size = base_size * 0.8, hjust=1, colour="grey50"))
  #������������Ϊ��׼��0.8������������Խڣ�x������ת90�ȣ�ɫ��Ϊ�л�



# ʹ��ggplot2��geom_tile����,����������ɫ����ͼ
# Ҳ���Ժ����ɵ�ʵ�ִ�ͳ�������ɫ����ƽ��䡣
p <- ggplot(data.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(low = "yellow",high = "red")

p + 
  theme_grey(base_size = base_size) + 
  labs(x = "",y = "") + 
  scale_x_continuous(expand = c(0, 0),labels=coln,breaks=1:length(coln)) + 
  scale_y_continuous(expand = c(0, 0),labels=rown,breaks=1:length(rown)) + 
  opts(axis.ticks = theme_blank(), 
       axis.text.x = theme_text(size = base_size * 0.8, angle = 90, hjust = 0, colour = "grey50"),
       axis.text.y = theme_text( size = base_size * 0.8, hjust=1, colour="grey50"))