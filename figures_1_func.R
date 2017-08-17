####################################### Functions from drawHm #############################
library(pheatmap)
library(gplots)

# 0.0 find the row index: rows with dist() val == NA
# input: the matrix
# output: matrix with those rows removed
matDelRows<-function(matr){
  distMt<-as.matrix(dist(matr))
  index<-which(is.na(distMt), arr.ind = T)
  index<-index[,1]
  index<-index[!duplicated(index)]
  matr<-matr[-index,]
  return(matr)
}

# 0.1 return the abs ceiling of floor of a digital numbel, depends on the sign
dpt<-function(x){
  if(x<0){
    x<-abs(floor(x))
  }
  else{
    x<-ceiling(x)
  }
  return(x)
}
# 0.2 Draw heatmap;
# mark=13:heatmap with unclustered rows and rownames;
# 23: unclustered rows and rownames and display numbers.
# 14, 24: unclusterd rows and cols
# 04: unclustered, no row name, no values
hm<-function(matrix, mark){
  depth<-max(dpt(min(matrix, na.rm = T)), dpt(max(matrix, na.rm = T)))
  c<-seq(-(depth+0.5),(depth), 0.5)
  len<-length(c)
  col<-bluered(len)
  # plotting
  if(mark==0) pheatmap(matrix, border_color = NA, breaks = c, color = col, show_rownames = F)
  else if(mark==1)   pheatmap(matrix, border_color = NA, breaks = c, color = col) # only row.names
  else if(mark==2) pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T) # TO SHOW THE NA
  else if(mark==13) pheatmap(matrix, border_color = NA, breaks = c, color = col, cluster_rows = F) # only row.names
  else if(mark==23) pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T, cluster_rows = F) # TO SHOW THE NA
  else if(mark==14) pheatmap(matrix, border_color = NA, breaks = c, color = col, cluster_rows = F, cluster_cols = F) # only row.names
  else if(mark==24) pheatmap(matrix, border_color = NA, breaks = c, color = col, display_numbers = T, cluster_rows = F, cluster_cols = F) # TO SHOW THE NA
  else if(mark==04) pheatmap(matrix, border_color = NA, breaks = c, color = col, show_rownames = F, cluster_rows = F, cluster_cols = F)
}

############################################################################################################
