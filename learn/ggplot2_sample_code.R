# ggplot2 codes from tutorial obtained by searching in the Internet

setwd("D:\\xnm\\work\\CODE\\learn")

 
# ��ͼ�Ĺؼ��Ǿ��࣬�������еķ����ǶԾ�������������ͽ����������ӻ���̶���
# ͨ�����plyr�������Ժܷ����ʵ�֡��������һ����Դ��WHO��������������ͼ�Ļ��ƽ��У�
# �������ݱ�׼������̬����Index��D��Ϊ�������˿����ݣ����������ٽ������ӻ���̶���
# ��geom_tile()������ͼ�Ļ��ƣ���ggplot2������ͨ��scale_fill_gradient2�����ֻ���ɫ���н��䡣
WHO<-read.csv("WHO.csv", header = TRUE)
require(plyr)
#�����˿�����������
WHO<-arrange(WHO, desc(D))
#�����ݵ�����ת��Ϊ���ӣ����̶����ĺõ�country��
#ͬ�����԰��վ���Ľ����������
WHO<- transform(WHO, Country = factor(Country, levels = unique(Country)))

require(reshape2)
require(ggplot2)
require(scales)
require(grid)
#melt����
m.WHO <- melt(WHO)
#��׼����ÿ������ӳ�䵽����Сֵ�����ֵӳ�䵽(0,1)����
m.WHO <- ddply(m.WHO, .(variable), transform, rescale = rescale(value))
#��׼������̬������
s.WHO <- ddply(m.WHO, .(variable), transform, rescale = scale(value))
require(ggplot2)
p<-ggplot(s.WHO, aes(variable, Country)) +
  #��tile�����л�����ͼ
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
    legend.key.size = unit(0.8, "cm"))#��Ҫ����grid��������legend�Ĵ�С



# ��ͼ��ɺ����һ������ͼƬ�������������ͼƬ�������������Ŀ��
# ������ȷ�������ʽ����ֱ�Ӳ��ý�ͼ�ķ�ʽ��ͼ���豸�н�ȡ���õ���ͼƬ�����ǵ��ӵġ�
# һ����������ͼƬӦ������ͼƬ�ߴ�������С������ʸ��ͼ���и�������Ⱦ������ν�Ŀ���ݡ�
# R����ͨ��֧��Cairoʸ��ͼ�δ�������⣬
# ���Դ�����������ʸ��ͼ��(PDF��PostScript��SVG) �� λͼ(PNG��JPEG�� TIFF)��
# ͬʱ֧���ں�̨�����и�������Ⱦ����ggplot2�ұȽ��Ƽ���ͼƬ�����ʽΪ����Cairo��������PDF��
# ��ΪPDF��ʽ���С��ͬʱ���Դ���Ϊ�����κθ�ʽ������ٽ�PDF����Ϊeps��ʽ����Photoshop�д������յĵ�����
# �������������ɫ�ʿռ��dpi��һ����־�ͳ�����Ҫ��dpi=300���ϣ��ȡ�
# ������Ҫע�����ggplot2�е������С���⣬��cookbook-rһ����ָ����
# ��ggplot2�о����������£�size�Ĵ�С��mm�ǣ���ϸ������Ҳ���Բο�stackover�����ۣ�
# ����theme()�ж�element_text()���size���е�������ʱ��size���԰�ֵ��points, pts�������б�ʾ��

# ������3��ggplot2�ֳ��õ�ͼƬ�����ʽ�����һ��������Ϊ20pts�������������Ϊ15pts��
# ��Ϊ80mm(3.15in)����Ϊ60mm(2.36in)��ͼΪ����
require(ggplot2)
require(Cairo)
ggplot() +
  geom_text(aes(x = 16, y = 16), label = "ABC", size = 11.28) + #�ߴ�Ϊ11.28mm����Ϊ32��
  geom_text(aes(x = 16, y = 14.5), label = "ABC", size = 32) + #�ߴ�Ϊ32mm
  labs( x = "x axis", y = "y axis") +
  ylim( c(14, 16.5)) +
  xlim( c(15.75, 16.25)) +
  theme(
    axis.title.x = element_text(size = 32),#�ߴ�Ϊ32��
    axis.title.y = element_text(size = 32))#�ߴ�Ϊ32��

x <- seq(-4,4, length.out = 1000)
y <-dnorm(x)
data <- data.frame(x, y)

#��Cairo�����
require(Cairo)
CairoPDF("plot1.pdf", 3.15, 3.15) #��λΪӢ��
ggplot(data, aes(x = x, y = y)) + geom_line(size = 1) +
  theme_bw()
dev.off() #�ر�ͼ���豸��ͬʱ����ͼƬ

plot2 <- ggplot(data, aes(x = x, y = y)) + geom_line(size = 1) +
  theme_bw()
#��ggsave�����Ĭ�ϼ�����Cairo�����п���ݴ���
ggsave("plot2.pdf", plot2, width = 3.15, height = 3.15) 

#RStudio���

# 
# ��extrafont���Ӣ������
# extrafont���ܹ�ֱ�ӵ��������ļ�����ͨ��Ghostscript(��Ҫ��װ����д�������������ɵ�PDF�У��������ɲο�������˵��
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
  geom_text(aes(x = 16, y = 15.50), label = "��������", size = 16,
            family = "STHUPO") +
  ylim(c(15.25, 16.50)) +
  labs(x = "", y = "") +
  theme_bw() #����RStudio���