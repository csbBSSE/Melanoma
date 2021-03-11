library(factoextra)
library(ggpubr)


#Cluster centres, for all genes 
mean_PCA <- function(data, cut, list){
  
  
  cluster <- cut
  axes = c(1,2)
  pca <- stats::prcomp(data[,list], scale = FALSE, center = FALSE)
  ind <- facto_summarize(pca, element = "ind", result = "coord", axes = 1:ncol(data))
  eig <- get_eigenvalue(pca)[axes,2]
  if(is.null(xlab)) xlab = paste0("Dim", axes[1], " (", round(eig[1],1), "%)") 
  if(is.null(ylab)) ylab = paste0("Dim", axes[2], " (", round(eig[2], 1),"%)")
  
  plot.data <- cbind.data.frame(ind, cluster = cluster, stringsAsFactors = TRUE)
  n <- ncol(plot.data)
  df <- plot.data[,-c(1,n-1,n)]
  
  
  mean <- data.frame()
  a <- unique(cluster)
  for (i in a) {
    
    mean[i,1:(n-3)] <- colMeans(df[plot.data$cluster==i,])
  }
  
  return(mean)
  
}



point_pairs_euc <- function(mean){
  
  mean <- as.data.frame(mean)
  n <- nrow(mean)
  pairs <- data.frame()
  l=1
  for (i in 1:(n-1)) {
    
    for (j in (i+1):n) {
      
      pairs[l,1:2] <- mean[i,1:2]
      pairs[l,3:4] <- mean[j,1:2]
      pairs[l,5] <- dist(mean[c(i,j),])
      pairs[l,6] <- i
      pairs[l,7] <- j
      l=l+1
    }
    
  }
  names(pairs) <- c("x","y","xend","yend", "dist", "From", "To")
  return(pairs)
  
}

avg <- function(x1,x2){
  
  x <- (x1+x2)/2
  return(x)
  
}

odd_even <- function(d){
  
  if(d%%2==0){
    
    return(-0.35)
    
  }else{
    
    return(+0.35)
    
  }
  
}

min_int_dist <- function(pairs){
  
  fr <- c(1,4)
  a <- data.frame()
  for (i in 1:2) {
    
    df <- pairs[(pairs$From== fr[i]|pairs$To== fr[i]),]
    a[i,1:4] <- df[which.min(df$dist),1:4]
    
  }
  return(a)
  
}

dist_ord <- function(d){
  
  df <- d[d$From== 1,]
  df <- dplyr::arrange(df, dist)
  df$To[df$To==2] <- c("T") 
  df$To[df$To==3] <- c("N") 
  df$To[df$To==4] <- c("U") 
  lab <- c(paste0("Distance from M: ", df$To[1]," < ",df$To[2]," < ",df$To[3] ))
  return(lab)
}

dist_PCA_plot <- function(data, cut, list){
  
  mean <- mean_PCA(data,cut$cluster)
  d <- as.data.frame(point_pairs_euc(mean))
  a <- PCA_all(data,list =  list,cut_f = cut)
  sol <- min_int_dist(d)
  
  for (i in 1:nrow(d)) {
    
    a <- a+ geom_segment(x=d$x[i],y=d$y[i],xend= d$xend[i], 
                         yend =d$yend[i],linetype =2, size=0.1)+
      geom_label(x= avg(d$xend[i],d$x[i]),  y= (avg(d$yend[i],d$y[i])+odd_even(i)), 
                 label =round(d$dist[i],2), label.size =0)
    
    
  }
  for (i in 1:nrow(sol)) {
    
    a <- a+ geom_segment(x=sol$x[i],y=sol$y[i],xend= sol$xend[i], 
                         yend =sol$yend[i], size=0.5)
    
    
  }
  print(a + geom_label(x= 0,  y= 3, label = dist_ord(d)))
  
  
}

clus_marker <- function(data, markers, seed =123, plot =T){
  
  set.seed(seed)
  cut <- kmeans(data[,markers],4)
  
  if(plot ==T){
    
    print(ComplexHeatmap::Heatmap(data[,markers], split = cut$cluster, column_order = markers,
                                  show_row_dend = F , heatmap_legend_param = list(title = "Scale"),
                                  show_row_names = F))
    
  }
  
  return(cut)
}

#PCA all genes
PCA_all <- function(data, markers, list, cut_f){
  
  require(factoextra)
  cut <- cut_f
  a <- fviz_cluster(cut, data = data[,list],xlab = "PC1", ylab = "PC2", geom = "point" ,
                    ellipse = F)+
    scale_shape_manual(values = c(16,16,16,16))+
    border()
  
  return(a)
  
}

list <- c("MITF", "FOS", "SMAD4", "STAT5A", "ETV5",  
          "SMAD3", "NR2F1", "NFIC",  "KLF4", "JUN","TFE3","NR3C1",
          "MAFB",  "TBX3","FOXF1","AHR", "TFAP2A")

data <- read.delim("Datasets/RACIPE.txt")
jpeg(paste0("Figures/Fig. 6/6A_Heatmap_clusters_RACIPE.jpeg"), width = 600, height = 520)
cut <- clus_marker(data, list,seed = 44)  
dev.off()

jpeg(paste0("Figures/Fig. 6/6A_trajectory_RACIPE.jpeg"), width = 560, height = 430)
dist_PCA_plot((data[,list]), cut, list)
dev.off()