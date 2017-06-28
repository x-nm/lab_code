
# ——————————————————————————————————————————————————————————————————————————————————————————————————
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256) # 深红深蓝，0不够明显不够白，
col = hmcol

col=heat.colors(256) # 不适合上下调的配色
col=topo.colors(100) # 颜色很丑
col=redgreen(75) # 适合上下调的色卡。但是有点丑

# ——————————————————————————————————————————————————————————————————————————————————————————————————
heatmap(mat1, na.rm = F, col = hmcol)


library("gplots")
# install.packages("gplots")
# heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,
#           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
heatmap.2(mat1, col=redgreen(75), scale="row",
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


# 想用不同的颜色来把样品组标记出来，那么我们可以使用ColSideColors参数来实现。
# 同时，我们希望变更热图的渐变填充色，可以使用col参数来实现。
color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
heatmap(data, col=topo.colors(100), ColSideColors=patientcolors, cexRow=0.5)





library(ggplot2)


##################################### test#######################
bmp(paste(title,".bmp",sep=""), width = wid, height = ht)
# plotting
heatmap.2(mat1, col = hmcol,scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none")

dev.off()
