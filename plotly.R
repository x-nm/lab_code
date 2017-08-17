# interactive ggplot
# 2017.07.11 by xnm

corrPlot<-function(dat,xcol,ycol){
  title <- paste0("Correlation of log2(norm+1) Expression of ",xcol,ycol)
  x = dat[,xcol]
  y = dat[,ycol]
  dat_new<-cbind(x,y,dat$ID)
  colnames(dat_new)<-c("x","y","ID")
  fit <- lm(y ~ x)
  R.sqr = summary(fit)$adj.r.squared
  xLoc <- min(x)+0.85*(max(x)-min(x))
  yLoc <- min(y)+0.2*(max(y)-min(y))
  p<-ggplot(as.data.frame(dat_new), aes(x = x, y = y, label = ID ))+
    geom_point(color = "slategray4",shape=1)+
    geom_smooth(method=lm,color="indianred4")+ # or indianred3
    labs(x=paste0(xcol,"_log2_norm_mean"),
         y=paste0(ycol,"_log2_norm_mean"))+
    ggtitle(title)+
    theme(
      plot.title = element_text(size = 1.2 * base_size, hjust = 0.5,face = "bold"),
      axis.text.x = element_text(vjust = 1, hjust = 1, size =  0.8 * base_size, color = "grey50"),
      axis.text.y = element_text(size = 0.8 * base_size, colour = "grey50")
    )+
    annotate("text", x=xLoc, y=yLoc, parse = T, label = paste("R^2 ==",round(R.sqr,4)))

  return(p)
  # ggsave(paste0(xcol,ycol,"_log2_norm.png"),width = 5, height = 4.8)
}

library(plotly)

p<-corrPlot(dat_merge,"W","O")

p<-corrPlot(WJHN,"W","H")
p<-ggplotly(p)
p

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
# chart_link = plotly_POST(p, filename="geom_point/stat-summary")
# chart_link



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
  dat_new<-as.data.frame(cbind(x,y,dat$ID))
  colnames(dat_new)<-c("x","y","ID")
  p<-ggplot(dat_new, aes(x = x, y = y, label=ID))+
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
  
  return(p)
}

p<-lfc_lgp(WH_deg,"WJ","DEG")
p<-ggplotly(p)
p
