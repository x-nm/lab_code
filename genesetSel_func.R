# genesetSel FUNCTION part
# for a collection of genesets, select genesets of interest
# 2017.6.7 by xnm


########################################################################################
# 0. small functions
########################################################################################

# 0.1
# Function: starLevel
# input: p value
# return: starLevel value
starLevel<-function(p){
  if(p>0.05) return("")
  else if(p<0.000001) return("****")
  else if(p<0.001) return("**")
  else return("*")
}


# 0.2
# Function: deg matched count for one sample
# input:  sample, dat_sub(geneset matched, fdr, fc selected)
# return: degCt for the sample
degMatchCt4sp<-function(sp, dat){
  return(nrow(dat[dat$sample == sp,]))
}

# 0.3
# Function: return a list of genes from a .gmx filename
# dependency: the global "genesedDir"
gainGeneset<-function(filename){
  genesetFile<-paste0(genesetDir,"\\",filename)
  gs<-readLines(genesetFile)
  gs<-gs[-c(1,2)]
  return(gs)
}


########################################################################################





########################################################################################
# 1. PCC
# 
# input: DEGct$degCount, degMatchCt for a geneset
# return: pcc, p, star
########################################################################################
pcc<-function(x,y){
 pcc<-cor(x,y,method = "pearson") 
 p<-cor.test(x,y,alternative = "two.sided",method = "pearson")$p.value
 if(is.na(p)) return(c(NA,NA,"")) # in case all samples matched 0
 star<-sapply(p, starLevel)
 return(c(pcc,p,star))
}
########################################################################################






########################################################################################
# 2. DISTRIBUTION TEST 
#
# input: matchedDEG count distribution of the geneset, prior distribution (prob)
# return: p-val, star-level(means H0 rejected or not, and reject level.)
########################################################################################
distriTest<-function(x,prior_prob){
  if(sum(x)==0) return(c(NA,""))
  
  chi<-chisq.test (x,
                   correct=TRUE,
                   p=prior_prob,
                   rescale.p= FALSE)
  p<-chi$p.value
  star<-sapply(p, starLevel)
  return(c(p, star))
}

########################################################################################







########################################################################################
# 3. ossification SCORE
#
# input: matchedDEG count, prior distribution (prob), index of interested samples
# return: star-type-score (in OB25, WHA, PO, anyone exceeds n*pi, score+=1)
########################################################################################
ossiScore<-function(matchedCt, prior_prob, index){
  score<-0
  n<-sum(matchedCt)
  for(i in 1:length(index)){
    id<-index[i]
    if(matchedCt[id]> n * prior_prob[id]) score<-score+1
  }
  star<-paste0(rep("*",score),collapse = "")
  return(star)
}
########################################################################################







# functions modified from newHeatmap_func.R
# 
# MAIN DIFFERENCE: 
#   "geneset" is a list of genes here, while it is a filename in the original function 
#
#   modified: dat_sub




########################################################################################
# 0.2 Whole Dataset preparation (with fdr, fc selected)
# return(dat_all)
########################################################################################
wholeDataSet<-function(fdr,fc,files){
  dataFile<-paste0("../fdr",fdr,"fc",fc,".txt")
  if(file.exists(dataFile)){
    dat_all<-read.csv(dataFile, header = T, sep = "\t")
  } 
  else{
    dat_all<-c()
    for(i in 1:length(files)){
      dat<-read.csv(paste0(files[i],"_LFC.txt"),header = T, sep = "\t")
      dat<-subset(dat,dat$padj<fdr)
      dat<-subset(dat,dat$LFC>log2(fc) | dat$LFC < log2((1/fc)))
      dat<-cbind(dat, rep(files[i],nrow(dat)))
      colnames(dat)[4]<-"sample"
      dat_all<-rbind(dat_all,dat)
    }
    write.table(dat_all,dataFile,quote = F, sep = "\t")
  }
  return(dat_all)
}
########################################################################################



########################################################################################
# 0.3 data subset
# 
# Dependency: global "samples"
# 
# input: gs means a list of genes
#
# return(dat_sub)
########################################################################################
# op=1, fix ID seq; (when geneset sequence matters more)
# op=2, fix gs seq; (when sample sequence matters more)

dataSubset<-function(dat, gs, op){
  
  dat<-dat[dat$ID%in%gs,] 
  
  # fix row/col sequence
  dat<-transform(dat, sample=factor(sample, levels = unique(samples)))
  if(op==1) return(transform(dat, ID=factor(ID, levels = unique(ID))))
  if(op==2) return(transform(dat, ID=factor(ID, levels = unique(gs))))
}
