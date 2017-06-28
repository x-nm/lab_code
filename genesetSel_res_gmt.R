# get res_all for .gmt


# gain geneset from .gmt
# "c5.all.v6.0.symbols" //5917, about 2min needed
# "c2.all.v6.0.symbols"
# "h.all.v6.0.symbols"
# one line for each geneset, \t, 1st name, 2nd discription(a url)

TIME1<-Sys.time()
res_all<-c()
genesetGroup<-c()
for(i in 1:nrow(gmt)){
  geneset<-as.vector(strsplit(gmt[i,1],split = "\t")[[1]])
  genesetGroup<-c(genesetGroup,geneset[1])
  
  geneset<-geneset[-c(1,2)] # geneset(a list) is the only key point for the consequent steps
  
  # sample distribution
  dat_sub<-dataSubset(dat_all, geneset, 1) # 1, sort as ID sequence; 2, sort as geneset sequence
  degMatchedCt<-sapply(samples,degMatchCt4sp,dat = dat_sub) # acturally, use as.vector(table(dat_sub$sample)) is enough
  # prior distribution
  distri_res<-distriTest(degMatchedCt,prior)
  # pcc
  pcc_res<-pcc(DEGct$degCount,degMatchedCt)
  # ossification score
  score<-ossiScore(degMatchedCt,prior, ossiIndex)
  # result
  res<-c(distri_res,score,pcc_res)
  res_all<-rbind(res_all,res)
}
colnames(res_all)<-c("distri_p","distri_rejection_level","ossi_score","pcc","pcc_p","pcc_level")
row.names(res_all)<-genesetGroup

TIME2<-Sys.time()
print(TIME2-TIME1)