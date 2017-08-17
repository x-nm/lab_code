#######################################################################################
# new heatmap for vc
#
# Functions part
#
# 1. separated data groups
# 2. match gene names (repliecated names allowed)
#       SPECIAL WAY TO SHOW THE REPLICATED ITEMS WITH DIFFERENT TIMES OF REP:
#       FIRST GET THE WHOLE SET: AABBCCD; THEN FILL IT BY SAMPLE, FILL IN BLANK IF NA.
# 3. no clustering of data and sample; show heatmap in the given order
#
# 
# INPUT:
# several group identifiers
# NORM. TXT: (Identifier) ID GROUP1 GROUP2 ...
# LFC.TXT: (Identifier) ID LFC padj
#
# * in LFC FILE, ID NA were excluded.
#
# working dir: where LFC.txt, samples.txt exists
# whole data file saved in ../ dir.
# geneset dir is defined in the main part DEF
# 
# 
# 2017.5.25 by xnm
# 
#######################################################################################

# setwd("D:\\xnm\\work\\17.5.25\\LFC/")

require(ggplot2)
require(grid)

#########################
# definition
#########################
base_size<-11
hColor<-"red"
mColor<-"white"
lColor<-"steelblue2"
#########################




########################################################################################
# 0. small functions:
########################################################################################

# 0.1.1 get geneset list
# input: geneset(FILENAME(no .gmx) for .gmx, or name of a geneset in .gmt file); type(gmt/gmx)
# dependency: global "genesetDir","gmt"data
# return: a list of genes, or geneset info
getGeneset<-function(geneset,type,rt="list"){
  if(type=="gmx"){
    genesetFile<-paste0(genesetDir,geneset,".gmx")
    gs<-readLines(genesetFile)
  }
  if(type=="gmt"){
    gs<-lineSplit(grep(geneset,gmt[,1],value = T))
  }
  info<-gs[c(1,2)]
  gs<-gs[-c(1,2)]
  if(rt=="list") return(gs)
  else if(rt=="info") return(info)
}


# 0.2 lineSplit, small one, in the aim of clear
# return: a vector of the line
lineSplit<-function(line, sep = "\t"){
  return(as.vector(strsplit(line,split = sep)[[1]]))
}


# 0.3
# Function: deg matched count for one sample
# input:  sample, dat_sub(geneset matched, fdr, fc selected)
# return: degCt for the sample
degMatchCt4sp<-function(sp, dat){
  return(nrow(dat[dat$sample == sp,]))
}

# 0.4
# Function: starLevel
# input: p value
# return: starLevel value
starLevel<-function(p){
  if(p>0.05 | p<=0) return("")
  else if(p<0.000001) return("****")
  else if(p<0.00001) return("***")
  else if(p<0.001) return("**")
  else return("*")
}

# 0.4.2
# star to number
star2num<-function(star){
  if(star=="") return(0)
  else if(star=="****") return(4)
  else if(star=="***") return(3)
  else if(star=="**") return(2)
  else return(1)
}

# 0.5 BOOL to star
bool2star<-function(bool){
  if(bool==TRUE) return("*")
  else return("")
}
########################################################################################




# borrowed from genesetSel_func.R

# 2. DISTRIBUTION TEST 
#
# input: matchedDEG count distribution of the geneset, prior distribution (prob)
# return: p-val, star-level(means H0 rejected or not, and reject level.)
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

# 3. ossification SCORE
#
# input: matchedDEG count, prior distribution (prob), index of interested samples
# return: star-type-score (in OB25, WHA, PO, anyone exceeds n*pi, score+=1)
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
# 1. heatmap with ggplot2
# 
# input: melted data, with col "LFC" for value, "sample" for groups, and "ID" for genes
#     (identifier)  ID  LFC  sample
# type: all (default), blank
# output: a ggplot2 object
# dependency: global "hgt", a dara.frame; global "geneset", global "samples"
########################################################################################
hm<-function(dat, type = "all", sqr = T, hgtStar = T){
  # level<-sapply(hgt[,4],star2num)
  # level<-as.data.frame(cbind(samples,level))
  # colnames(level)<-c("sample","starLevel")
  if(hgtStar == T){
    spName<-paste(hgt[,4],samples)
  }
  else{
    spName<-samples
  }
  
  
  # ggplot2
  plot<-ggplot(dat, aes(sample, ID, label = dup)) +
    geom_tile(aes(fill=LFC)) +
    scale_fill_gradient2(mid=mColor, high=hColor, low=lColor, name = "LogFC") +
    # labs(x="Sample", y="Gene", face = "bold") +
    labs(x="",y="")+
    theme_bw() +
    theme(
      plot.title = element_text(size = 1.2*base_size, hjust = 0.5,face = "bold"),
      # axis.title.x=element_text(size = base_size),
      # axis.title.y = element_text(size = base_size),
      axis.text.x = element_text(angle = 90,  vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
      axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50"),
      legend.title = element_text(size = 0.8 *  base_size),
      legend.text = element_text(size = 0.8 *  base_size),
      legend.key.size = unit(0.8, "cm")
      # axis.line.x = element_line(color = level$starLevel,linetype = 2, size = 0.2 * base_size)
      )+
    geom_text(size = 0.3 * base_size, colour = "grey50")+
    ggtitle(geneset)+ # title
    scale_x_discrete(breaks=samples, labels=spName) # x label:spname+stars
  
  if(type=="blank"){
    plot<-plot+
      ylab("")+
      theme(axis.text.y = element_blank())
  }
  if(sqr==T){
    plot<-plot+coord_fixed() # square box
  }
  
  return(plot)
}
########################################################################################



########################################################################################
# 2. Whole Dataset preparation (with fdr, fc selected)
# 
# if aggr = T, aggregate duplicated LFCs (after DEG selection)
# return(dat_all)
########################################################################################
wholeDataSet<-function(fdr,fc,files, aggr = T){
  if(aggr){
    dataFile<-paste0("../fdr",fdr,"fc",fc,"_dup_aggr.txt")
  }
  else{
    dataFile<-paste0("../fdr",fdr,"fc",fc,".txt")
  }
  
  if(file.exists(dataFile)){
    dat_all<-read.csv(dataFile, header = T, sep = "\t")
  } 
  else{
    dat_all<-c()
    for(i in 1:length(files)){
      dat<-read.csv(paste0(files[i],"_LFC.txt"),header = T, sep = "\t")
      dat<-subset(dat,dat$padj<fdr)
      # dat<-subset(dat,dat$LFC>log2(fc) | dat$LFC < log2((1/fc)))
      dat<-cbind(dat, rep(files[i],nrow(dat)))
      colnames(dat)[4]<-"sample"
      
      if(aggr){
        dat<-cbind(dat,rep(F,nrow(dat)))
        colnames(dat)[5]<-"dup"
        dupIndex<-which(duplicated(dat$ID))
        
        if(length(dupIndex)!=0){
          dupDat<-dat[dupIndex,]
          dat<-dat[-dupIndex,]
          for(j in 1:nrow(dupDat)){
            dat[which(dat$ID==dupDat$ID[j]),"LFC"]<-dat[which(dat$ID==dupDat$ID[j]),"LFC"][1]+dupDat$LFC[j]
            dat[which(dat$ID==dupDat$ID[j]),"dup"]<-T
          }
        }
        
        dat$dup<-sapply(dat$dup,bool2star)
      }
      dat<-subset(dat,dat$LFC>log2(fc) | dat$LFC < log2((1/fc)))
      dat_all<-rbind(dat_all,dat)
    }
    write.table(dat_all,dataFile,quote = F, sep = "\t")
  }
  return(dat_all)
}


# ori_whole_data_dup_del
origWhole<-function(files){
  dat_all<-c()
  for(i in 1:length(files)){
    dat<-read.csv(paste0(files[i],"_LFC.txt"),header = T, sep = "\t")
    dat<-cbind(dat, rep(files[i],nrow(dat)))
    colnames(dat)[4]<-"sample"
    
    dupIndex<-which(duplicated(dat$ID))
    if(length(dupIndex)!=0){
      dat<-dat[-dupIndex,]
    }
    dat_all<-rbind(dat_all,dat)
  }
  write.table(dat_all,"ori_whole_data_dup_del",quote = F, sep = "\t")
}

########################################################################################



########################################################################################
# 3. data subset
# 
# Dependency: global "samples",global "gs"(for function2)
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

dataSubset4sps<-function(dat, sps,op=1){
  dat<-dat[dat$sample%in%sps,] 
  # fix row/col sequence
  dat<-transform(dat,sample=factor(sample, levels = unique(sps)))
  if(op==1) return(transform(dat, ID=factor(ID, levels = unique(ID))))
  if(op==2) return(transform(dat, ID=factor(ID, levels = unique(gs))))
}
########################################################################################




########################################################################################
# 4. save 
# 
# Dependency: global variable: dat_sub
########################################################################################
hmSave<-function(plot, title, a = 4.2, b = 3, save.ht = 5.5){
  # save.ht 3.15~5.5; a 4.2, b 3
  n<-length(levels(dat_sub$ID))
  m<-length(levels(dat_sub$sample))
  a<- a/(n+3) * save.ht # a=4.2x, x=h/(n+3); 
  b<- b/(n+3) * save.ht # b=3x
  save.wid<-m/n * (save.ht - b) + 2*a
  ggsave(paste0("output/",title), plot, width = save.wid, height = save.ht) 
}
########################################################################################




########################################################################################
# 5. hypergeometric distribution test, for each sample
#
# input: 
#    melt WHOLE data of the sample, (ID  LFC padj  sample)
#    geneset (list)
#    fdr, fc (for find the DEGs)
# output:
#    p-val       
    
########################################################################################
hyperGeomTest<-function(dat, gs, fdr, fc){
  gsSize<-length(gs[gs%in%dat$ID]) # meaningful geneset: with GENE matches the dataset$ID
  
  bgSize<-nrow(dat)
  
  deg<-subset(dat,dat$padj<fdr)
  deg<-subset(deg,deg$LFC>log2(fc) | deg$LFC < log2((1/fc)))
  DEGsize<-nrow(deg)
  
  overlap<-length(deg[deg$ID%in%gs,"ID"])
  
  pval<-phyper(q = overlap, 
               m = DEGsize , # number of red balls
               n = bgSize, # number of other balls
               k = gsSize, # number of balls drawn
               lower.tail = F) # hypergeometric test, SIGMA x>q
  return(pval)
}
########################################################################################



########################################################################################
# 6. hyperGeomTest for all samples
# 
# input: original whole data, samples (a vector of samples used, in the fixed order), geneset(list), fdr, fc
# output: a vector of hyperGeom p-vals, in the order of samples
########################################################################################

HGTpval4samples<-function(dat, sampleList, geneset, fdr, fc){
  pval<-c() # vector for return
  for(i in 1:length(sampleList)){
    spData<-dat[dat$sample==sampleList[i],]
    p<-hyperGeomTest(spData, geneset, fdr, fc)
    pval<-c(pval,p)
  }
  
  padj<-p.adjust(pval,"fdr",n = 8000) # ADJUST; add an n=6000 to restrict the hgt levels
  star<-sapply(padj,starLevel)
  
  pval<-cbind(sampleList,pval,padj, star)
  colnames(pval)<-c("samples","pval","fdr","star")
  return(pval)
}
########################################################################################




########################################################################################
# 7. Distribution plot
# 
# "prob density" plot for prior distribution and post distribution in one plot
# 
# input: prior ditribution(prob vector) , post distribution COUNT(NAMED prob vector)
# depencecy: global "geneset", "samples","ossiIndex",  global "hgt" 
#     function distriTest, ossiScore, starLevel..
########################################################################################


distriPlot<-function(prior, post_ct){
  post<-post_ct/sum(post_ct)
  dat<-c(prior,post)
  n<-length(prior)
  cls<-c(rep("prior",n),rep("post",n))
  sps<-c(names(post),names(post))
  dat<-cbind(dat,sps,cls)
  colnames(dat)<-c("prob","samples","class")
  dat<-as.data.frame(dat)
  # fix with samples sequence
  dat<-transform(dat, samples=factor(samples, levels = unique(samples)))

  dis_p<-distriTest(post_ct,prior)
  dis_p<-paste("p = ", dis_p[2])
  score<-ossiScore(post_ct,prior,ossiIndex)
  score<-paste("OssiScore = ", score)
  
  spName<-paste(hgt[,4],samples) # change x-axis label with hgt p-val
  
  p<-ggplot(data=dat, aes(x=samples,y=prob,color=class,group=class,shape = class)) +
    ylab("Probability")+
    xlab("")+
    geom_line() +
    geom_point()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,  vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 1.2*base_size, face = "bold", hjust = 0.5)
    )+
    annotate("text",x = 15 , y = 25, label = score) +
    annotate("text", x = 15, y = 23, label = dis_p) +
    ggtitle(geneset)+ # title
    scale_x_discrete(breaks=samples, labels=spName) # x label:spname+stars
  
  # p
  
  return(p)
}


# post_ct<-degMatchedCt
# prior<-prior_prob


########################################################################################













