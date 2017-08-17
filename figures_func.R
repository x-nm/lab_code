# functions for drawing small figures in rabbit project
# 2017.07.10 by xnm


##########################
base_size<-11
hColor<-"red"
mColor<-"white"
lColor<-"steelblue2"
##########################


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                           FUNCTIONS FOR 0.2                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(ggplot2)

# 2.1 aggrDup
# input: norm dat with an ID col
aggrDup<-function(dat){
  dupIndex<-which(duplicated(dat$ID))
  if(length(dupIndex)!=0){
    dupDat<-dat[dupIndex,]
    dat<-dat[-dupIndex,]
    for(j in 1:nrow(dupDat)){
      dat[which(dat$ID==dupDat$ID[j]),-1]<-dat[which(dat$ID==dupDat$ID[j]),-1]+dupDat[j,-1]
    }
  }
  return(dat)
}

# 2.2 normMean
normMean<-function(dat,index1,index2,log2=T, aggr=F){
  if(aggr){
    dat<-aggrDup(dat)
  }
  num_dat<-dat[,index1:index2]
  if(log2){
    num_dat<-log2(num_dat+1)
  }
  res<-as.data.frame(cbind(dat[,1],apply(num_dat,1,mean)))
  colnames(res)<-c("ID","Mean_Of_Norm")
  return(res)
}

# 2.3.1 ID
id<-function(id,bool){
  if(bool){
    return(id)
  }
  else{
    return(NA)
  }
}

# 2.3 corrPlot
# input:
#   dat: ID W J H N O
#   xcol: "W"
#   ycol: "H"
corrPlot<-function(dat,xcol,ycol,type="norm"){
  x = dat[,xcol]
  y = dat[,ycol]
  # index<-dat[,xcol]>15 & dat[,ycol]>15
  # ID<-sapply(dat$ID,id,bool=index)
  # n<-nrow(dat)
  # ID<-c()
  # for(i in 1:n){
    # ID[i]<-id(dat$ID[i],index[i])
  # }
  # dat_new<-as.data.frame(cbind(x,y,ID))
  # colnames(dat_new)<-c("x","y","ID")
  dat_new<-as.data.frame(cbind(x,y))
  colnames(dat_new)<-c("x","y")
  fit <- lm(y ~ x,na.action = na.exclude)
  R.sqr = summary(fit)$adj.r.squared
  xLoc <- min(x)+0.85*(max(x)-min(x))
  yLoc <- min(y)+0.2*(max(y)-min(y))

  if(type=="norm"){
    title <- paste0("Correlation of log2(norm+1) Expression of ",xcol,ycol)
    ggplot(dat_new, aes(x = x, y = y))+
      geom_point(color = "slategray4",shape=1)+ # or no shape = 1
      geom_smooth(method=lm,color="indianred4")+ # or indianred3
      labs(x=paste0(xcol,"_log2_norm_mean"),
           y=paste0(ycol,"_log2_norm_mean"))+
      ggtitle(title)+
      theme(
        plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
        axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
      )+
      # geom_text(size = 0.5 * base_size, colour = "grey50", position = "identity", na.rm = T)+
      annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))
    ggsave(paste0(xcol,ycol,"_log2_norm.png"),width = 5, height = 4.8)
  }
  else if(type=="lfc"){
    title <- paste0("Correlation of LogFC of ",xcol,ycol)
    ggplot(dat_new, aes(x = x, y = y))+
      geom_point(color = "slategray4")+ # or no shape = 1
      geom_smooth(method=lm,color="indianred4")+ # or indianred3
      labs(x=paste0(xcol,"_lfc"),
           y=paste0(ycol,"_lfc"))+
      ggtitle(title)+
      theme(
        plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
        axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
      )+
      # geom_text(size = 0.5 * base_size, colour = "grey50", position = "identity", na.rm = T)+
      annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))
    ggsave(paste0(xcol,ycol,"_lfc_corr.png"),width = 5, height = 4.8)
  }else if(type=="lfc_deg"){
    title <- paste0("Correlation of LogFC of DEGs of ",xcol,ycol)
    ggplot(dat_new, aes(x = x, y = y))+
      geom_point(color = "slategray4")+ # or no shape = 1
      geom_smooth(method=lm,color="indianred4")+ # or indianred3
      labs(x=paste0(xcol,"_lfc"),
           y=paste0(ycol,"_lfc"))+
      ggtitle(title)+
      theme(
        plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
        axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
      )+
      # geom_text(size = 0.5 * base_size, colour = "grey50", position = "identity", na.rm = T)+
      annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))
    ggsave(paste0(xcol,ycol,"_lfc_deg_corr.png"),width = 5, height = 4.8)
  }
  
}




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                           FUNCTIONS FOR 0.3                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
dat_DEG<-function(dat, fdr, fc){
  dat<-subset(dat,dat$padj<fdr)
  dat<-subset(dat,dat$LFC>log2(fc) | dat$LFC < log2((1/fc)))
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                           FUNCTIONS FOR 0.4                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# (WH_deg, "WJ", "DEG")
lfc_lgp<-function(dat,xcol,type="all"){
  if(type=="all"){
    plot_title<-paste0("LFC-log(p) plot of ",xcol)
    save_title<-paste0(xcol,"_lfc_lgp.png")
  }else if(type=="DEG"){
    plot_title<-paste0("LFC-log(p) plot of DEG of ",xcol)
    save_title<-paste0(xcol,"_deg_lfc_lgp.png")
  }
  ycol<-paste0(xcol,".p-val")
  x = dat[,xcol]
  y = dat[,ycol]
  y = -log10(y)
  dat_new<-as.data.frame(cbind(x,y))
  colnames(dat_new)<-c("x","y")
  ggplot(dat_new, aes(x = x, y = y))+
    geom_point(color = "slategray4")+ # or no shape = 1
    geom_abline(slope = 0,intercept = 1.3, color = "indianred4")+
    labs(x=paste0(xcol,"_lfc"),
         y=paste0(xcol,"_p-val"))+
    ggtitle(plot_title)+
    theme(
      plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
      axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
      axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
    )
  
  ggsave(save_title,width = 5, height = 4.8)
  
  
  
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                           FUNCTIONS FOR 0.5                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(VennDiagram)


# 2. Whole Dataset preparation (with fdr, fc selected) 
# if aggr = T, aggregate duplicated LFCs (after DEG selection)
# return(dat_all)
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
      dat<-subset(dat,dat$LFC>log2(fc) | dat$LFC < log2((1/fc)))
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
      dat_all<-rbind(dat_all,dat)
    }
    write.table(dat_all,dataFile,quote = F, sep = "\t")
  }
  return(dat_all)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                           FUNCTIONS FOR TOP25                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# dat: ID, Value, Sample
simple_heatmap<-function(dat,title,vec = T){
  if(vec){
    p<-ggplot(dat,aes(x=Sample,y=ID))
  }else{
    p<-ggplot(dat,aes(x=ID,y=Sample))
  }
  p<-p+
    geom_tile(aes(fill=Value))+
    scale_fill_gradient2(mid=mColor, high=hColor, low=lColor) +
    xlab("")+
    ylab("")+
    theme_bw() +
    theme(
      plot.title = element_text(size = 1.2*base_size, hjust = 0.5,face = "bold"),
      axis.text.x = element_text(angle = 90,  vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
      axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50"),
      legend.title = element_text(size = 0.8 *  base_size),
      legend.text = element_text(size = 0.8 *  base_size),
      legend.key.size = unit(0.8, "cm")
    )+
    ggtitle(title)+ # title
    coord_fixed()
  return(p)
}


# Dependency: global "saveDir"
top_norm<-function(dat,sp,savelb,thres=25){
  dat<-dat[order(dat[,sp],decreasing = T),]
  dat<-dat[1:25,c("ID",sp)]
  colnames(dat)<-c("ID","Norm.Mean")
  saveTitle<-paste0(saveDir,"top25_",savelb,"_",sp,".txt")
  write.table(dat,saveTitle,quote = F, sep="\t",row.names = F)
  return(dat)
}

# Dependency: global "saveDir"
top_deg<-function(dat,up,sp,savelb,thres=25){
  if(up){
    dat<-dat[order(dat[,"LFC"],decreasing = T),]
  }else{
    dat<-dat[order(dat[,"LFC"],decreasing = F),]
  }
  dat<-dat[1:25,c("ID","LFC")]
  saveTitle<-paste0(saveDir,"top25_",savelb,"_",sp,".txt")
  write.table(dat,saveTitle,quote = F, sep="\t",row.names = F)
  return(dat)
}

# Dependency: global "saveDir"
top_deg_pval<-function(dat,up,sp,savelb,thres=25){
  if(up){
    dat<-subset(dat,dat$LFC>0)
    dat<-dat[order(dat[,"padj"],decreasing = F),]
  }else{
    dat<-subset(dat,dat$LFC<0)
    dat<-dat[order(dat[,"padj"],decreasing = F),]
  }
  dat<-dat[1:25,c("ID","LFC")]
  saveTitle<-paste0(saveDir,"top25_",savelb,"_",sp,".txt")
  write.table(dat,saveTitle,quote = F, sep="\t",row.names = F)
  return(dat)
}

