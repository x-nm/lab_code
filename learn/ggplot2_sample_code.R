# ggplot2 codes from tutorial obtained by searching in the Internet

setwd("D:\\xnm\\work\\CODE\\learn")

 
# 热图的关键是聚类，两个可行的方案是对聚类结果进行排序和将聚类结果因子化后固定，
# 通过结合plyr包，可以很方便的实现。这里采用一组来源于WHO国家数据来对热图的绘制进行，
# 首先数据标准化和正态化后按Index的D（为各国的人口数据）进行排序，再将其因子化后固定，
# 用geom_tile()进行热图的绘制，在ggplot2种已能通过scale_fill_gradient2在三种基本色进行渐变。
WHO<-read.csv("WHO.csv", header = TRUE)
require(plyr)
#按总人口数排列数据
WHO<-arrange(WHO, desc(D))
#将数据的名字转换为因子，并固定已拍好的country，
#同理可以按照聚类的结果进行排列
WHO<- transform(WHO, Country = factor(Country, levels = unique(Country)))

require(reshape2)
require(ggplot2)
require(scales)
require(grid)
#melt数据
m.WHO <- melt(WHO)
#标准化，每排数据映射到按最小值和最大值映射到(0,1)区间
m.WHO <- ddply(m.WHO, .(variable), transform, rescale = rescale(value))
#标准化并正态化数据
s.WHO <- ddply(m.WHO, .(variable), transform, rescale = scale(value))
require(ggplot2)
p<-ggplot(s.WHO, aes(variable, Country)) +
  #用tile来进行绘热力图
  geom_tile(aes(fill=rescale)) +
  scale_fill_gradient2(mid="black", high="red", low="green", name = "Intensity") +
  labs(x="Country", y="Index", face = "bold") +
  theme_bw() +
  theme(
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=12, colour="grey50"),
    axis.text.y=element_text(size=12, colour="grey50"),
    legend.title=element_text(size=14),
    legend.text=element_text(size=12),
    legend.key.size = unit(0.8, "cm"))#需要载入grid包来调整legend的大小



# 绘图完成后最后一步便是图片输出，高质量的图片输出让人赏心悦目，
# 而不正确的输出方式或者直接采用截图的方式从图形设备中截取，得到的图片往往是低劣的。
# 一幅高质量的图片应当控制图片尺寸和字体大小，并对矢量图进行高质量渲染，即所谓的抗锯齿。
# R语言通过支持Cairo矢量图形处理的类库，
# 可以创建高质量的矢量图形(PDF，PostScript，SVG) 和 位图(PNG，JPEG， TIFF)，
# 同时支持在后台程序中高质量渲染。在ggplot2我比较推荐的图片输出格式为经过Cairo包处理的PDF，
# 因为PDF格式体积小，同时可以储存为其他任何格式，随后再将PDF储存为eps格式并在Photoshop中打开做最终的调整，
# 例如调整比例、色彩空间和dpi（一般杂志和出版社要求dpi=300以上）等。
# 额外需要注意的是ggplot2中的字体大小问题，在cookbook-r一书中指出，
# 在ggplot2中绝大多数情况下，size的大小以mm记，详细的讨论也可以参考stackover的讨论，
# 而在theme()中对element_text()里的size进行调整，此时的size是以磅值（points, pts）来进行表示。

# 下面以3种ggplot2种常用的图片输出方式，输出一幅主标题为20pts，横纵坐标标题为15pts，
# 长为80mm(3.15in)，宽为60mm(2.36in)的图为例。
require(ggplot2)
require(Cairo)
ggplot() +
  geom_text(aes(x = 16, y = 16), label = "ABC", size = 11.28) + #尺寸为11.28mm，即为32磅
  geom_text(aes(x = 16, y = 14.5), label = "ABC", size = 32) + #尺寸为32mm
  labs( x = "x axis", y = "y axis") +
  ylim( c(14, 16.5)) +
  xlim( c(15.75, 16.25)) +
  theme(
    axis.title.x = element_text(size = 32),#尺寸为32磅
    axis.title.y = element_text(size = 32))#尺寸为32磅

x <- seq(-4,4, length.out = 1000)
y <-dnorm(x)
data <- data.frame(x, y)

#用Cairo包输出
require(Cairo)
CairoPDF("plot1.pdf", 3.15, 3.15) #单位为英寸
ggplot(data, aes(x = x, y = y)) + geom_line(size = 1) +
  theme_bw()
dev.off() #关闭图像设备，同时储存图片

plot2 <- ggplot(data, aes(x = x, y = y)) + geom_line(size = 1) +
  theme_bw()
#用ggsave输出，默认即以用Cairo包进行抗锯齿处理
ggsave("plot2.pdf", plot2, width = 3.15, height = 3.15) 

#RStudio输出

# 
# 用extrafont输出英文字体
# extrafont包能够直接调用字体文件，再通过Ghostscript(需要安装）将写入的字体插入生成的PDF中，具体代码可参考了作者说明
# https://cos.name/2014/01/showtext-interesting-fonts-and-graphs/

#showtext
require(showtext)
require(ggplot2)
require(Cairo)
font.add("BlackoakStd", "C://Windows//Fonts//BlackoakStd.otf")
font.add("BrushScriptStd", "C://Windows//Fonts//BrushScriptStd.otf")
font.add("times", "C://Windows//Fonts//times.ttf")
font.add("STHUPO", "C://Windows//Fonts//STHUPO.ttf")
CairoPDF("showtext_output", 8, 8)
showtext.begin()
ggplot() +
  geom_text(aes(x = 16, y = 16.25), label = "Blackoak Std", size = 8, 
            family = "BlackoakStd") +
  geom_text(aes(x = 16, y = 16), label ="Brush Script Std", size = 16,
            family = "BrushScriptStd") +
  geom_text(aes(x = 16, y = 15.75), label = "Times New Roman", size = 16,
            family = "times") +
  geom_text(aes(x = 16, y = 15.50), label = "华文琥珀", size = 16,
            family = "STHUPO") +
  ylim(c(15.25, 16.50)) +
  labs(x = "", y = "") +
  theme_bw() #在用RStudio输出