# get the result for .gmx in a dir
res_all<-c()
for(i in 1:length(genesetGroup)){
  geneset<-gainGeneset(genesetGroup[i]) # geneset(a list) is the only key point for the consequent steps
  
  # sample distribution
  dat_sub<-dataSubset(dat_all, geneset, 1) # 1, sort as ID sequence; 2, sort as geneset sequence
  degMatchedCt<-sapply(samples,degMatchCt4sp,dat = dat_sub) # acturally, use as.vector(table(dat_sub$sample)) is enough
  # distribution test
  distri_res<-distriTest(degMatchedCt,prior)
  # pcc
  pcc_res<-pcc(DEGct$degCount,degMatchedCt)
  # ossification score
  ossiIndex<-match(ossiSp,samples)
  score<-ossiScore(degMatchedCt,prior, ossiIndex)
  # result
  res<-c(distri_res,score,pcc_res)
  res_all<-rbind(res_all,res)
}
colnames(res_all)<-c("distri_p","distri_rejection_level","ossi_score","pcc","pcc_p","pcc_level")
row.names(res_all)<-genesetGroup