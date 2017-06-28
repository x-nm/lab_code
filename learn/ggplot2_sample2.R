library(ggplot2)
hc<-hclust(dist(data))
rowInd<-hc$order
hc<-hclust(dist(t(data)))
colInd<-hc$order
data.m<-data[rowInd,colInd] #聚类分析的作用是为了色块集中，显示效果好。如果本身就对样品有分组，
# 基因有排序，就可以跳过这一步。

data.m<-apply(data.m,1,rescale) #以行为基准对数据进行变换，使每一行都变成［0,1］之间的数字。
# 变换的方法可以是scale,rescale等等，按照自己的需要来变换。

data.m<-t(data.m) #变换以后转置了。
coln<-colnames(data.m) 
rown<-rownames(data.m) #保存样品及基因名称。因为geom_tile会对它们按坐标重排，
# 所以需要使用数字把它们的序列固定下来。

colnames(data.m)<-1:ncol(data.m)
rownames(data.m)<-1:nrow(data.m)
data.m<-melt(data.m) #转换数据成适合geom_tile使用的形式
head(data.m)
# X1 X2     value
# 1  1  1 0.1898007
# 2  2  1 0.6627467
# 3  3  1 0.5417057
# 4  4  1 0.4877054
# 5  5  1 0.5096474
# 6  6  1 0.2626248

base_size<-12 #设置默认字体大小，依照样品或者基因的多少而微变。


p <- ggplot(data.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  #设定横坐标为以前的列，纵坐标为以前的行，填充色为转换后的数据
  scale_fill_gradient(low = "white",high = "steelblue")
  #设定渐变色的低值为白色，变值为钢蓝色。


p + 
  theme_grey(base_size = base_size) + 
  labs(x = "", #设置xlabel及ylabel为空
       y = "") + 
  scale_x_continuous(expand = c(0, 0),labels=coln,breaks=1:length(coln)) + 
  #设置x坐标扩展部分为0，刻度为之前的样品名
  scale_y_continuous(expand = c(0, 0),labels=rown,breaks=1:length(rown)) + 
  #设置y坐标扩展部分为0，刻度为之前的基因名
  opts(axis.ticks = theme_blank(), 
       axis.text.x = theme_text(size = base_size * 0.8, angle = 90, hjust = 0, colour = "grey50"), 
       axis.text.y = theme_text(size = base_size * 0.8, hjust=1, colour="grey50"))
  #设置坐标字体为基准的0.8倍，贴近坐标对节，x坐标旋转90度，色彩为中灰



# 使用ggplot2中geom_tile函数,钢蓝渐白配色的热图
# 也可以很轻松的实现传统渐变填充色，红黄渐变。
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