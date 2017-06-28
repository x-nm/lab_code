# get all the descriptions of .gmx
# 2017.6.8 by xnm

genesetDir<-"D:\\xnm\\work\\17.2.15\\KeyGene\\geneset"

gainGsDes<-function(filename){
  genesetFile<-paste0(genesetDir,"\\",filename)
  gs<-readLines(genesetFile)
  gs<-gs[c(1,2)]
  return(gs)
}

genesetGroup<-grep(".gmx$",list.files(genesetDir),value = T)

des<-c("geneset","description")
for(i in 1:length(genesetGroup)){
  geneset<-gainGsDes(genesetGroup[i])
  des<-rbind(des, geneset)
}
row.names(des)<-c("",genesetGroup)

write.table(des,"des_gmx.txt",quote = F, sep = "\t")
