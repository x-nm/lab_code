
# ����������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256) # ���������0�������Բ����ף�
col = hmcol

col=heat.colors(256) # ���ʺ����µ�����ɫ
col=topo.colors(100) # ��ɫ�ܳ�
col=redgreen(75) # �ʺ����µ���ɫ���������е��

# ����������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������
heatmap(mat1, na.rm = F, col = hmcol)


library("gplots")
# install.packages("gplots")
# heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,
#           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
heatmap.2(mat1, col=redgreen(75), scale="row",
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


# ���ò�ͬ����ɫ������Ʒ���ǳ�������ô���ǿ���ʹ��ColSideColors������ʵ�֡�
# ͬʱ������ϣ�������ͼ�Ľ������ɫ������ʹ��col������ʵ�֡�
color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
heatmap(data, col=topo.colors(100), ColSideColors=patientcolors, cexRow=0.5)





library(ggplot2)


##################################### test#######################
bmp(paste(title,".bmp",sep=""), width = wid, height = ht)
# plotting
heatmap.2(mat1, col = hmcol,scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none")

dev.off()