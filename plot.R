# plot other plots needed
# 2017.2.21 by xnm

setwd("D:\\xnm\\work\\17.2.15\\KeyGene")

library(pheatmap)


############################################## CODE BANK############################

# ������������ͼ���Զ���ͼƬ����С��΢�е㲻̫��
pheatmap(mat, border_color = NA, show_rownames = F ,filename = "test.png")


# �ԳƵ�ɫ�������֮ǰ��pheatmap��bmpҲ�ܺ���; ��Ϊ��scale="row"���������Եø����ԣ�
library("gplots")
heatmap.2(mat1, col=redgreen(75), scale="row",
          key=TRUE, symkey=FALSE, density.info="none", trace="none")

# ԭ������kegGene.R���pheatmap
pheatmap(mat1, border_color = NA, display_numbers = T)
# scale��pheatmap
pheatmap(mat1, border_color = NA, scale = "row") # ��Ȼ���öԱ��Ը�������һ�㡣�������ɫ��0��ø����Կ��ܸ��ã�
# ����scale֮������µ�ֵҲ��仯��ֻ����Ե���-1~1�ķ�Χ�ڶ��ѡ��������ޡ�

# ����breaks���heatmap(0Ϊ��ɫ)
library(gplots)
depth<-round(max(abs(min(mat1)),abs(max(mat1))))
c<-c(-(depth+1):depth)
len<-length(c)
# a<-colorpanel(len,low="blue",mid = "white", high = "red")
b<-bluered(len)
pheatmap(mat1, border_color = NA, breaks = c, color = b)
pheatmap(mat1, border_color = NA, breaks = c, color = b, display_numbers = T) # TO SHOW THE NA



####################################################################################

# 1. ehm of all data



# ehm of one sample's DEG; as cannot allocate data with too many items in this PC
mat<-read.table("PO_NORM.txt", header = T, row.names = 1)
ind<-read.table("PO_LFC.txt",header = T)
ind<-subset(ind, ind$padj<0.05)
index<-na.omit(match(row.names(mat),ind$ID))
index<-match(ind$ID,row.names(mat))
mat<-mat[index,]

# log2(x)
mat2<-log2(mat)
mat2[mat2==-Inf]<-0
mat2[mat2==0]<-min(mat2)-1
sum(is.na(mat2))
pheatmap(mat2,fontsize = 5,show_rownames = F)

# log2(x+1)
mat<-log2(mat+1)
mat<-as.matrix(mat)
sum(is.na(mat))
pheatmap(mat,fontsize = 5,show_rownames = F)


# 2. lhm of all data
mat<-read.table("WHPA_LFC.txt",header = T)



# 3. dlhm of all data