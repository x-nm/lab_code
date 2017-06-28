# plot other plots needed
# 2017.2.21 by xnm

setwd("D:\\xnm\\work\\17.2.15\\KeyGene")

library(pheatmap)


############################################## CODE BANK############################

# 画表达量的热图，自动存图片，大小稍微有点不太好
pheatmap(mat, border_color = NA, show_rownames = F ,filename = "test.png")


# 对称的色卡，替代之前的pheatmap画bmp也很合适; 因为有scale="row"所以容易显得更明显？
library("gplots")
heatmap.2(mat1, col=redgreen(75), scale="row",
          key=TRUE, symkey=FALSE, density.info="none", trace="none")

# 原来用在kegGene.R里的pheatmap
pheatmap(mat1, border_color = NA, display_numbers = T)
# scale的pheatmap
pheatmap(mat1, border_color = NA, scale = "row") # 果然看得对比性更加明显一点。如果换个色盘0变得更明显可能更好？
# 但是scale之后的上下调值也会变化！只是相对到了-1~1的范围内而已。不可以噢。

# 加入breaks后的heatmap(0为白色)
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
