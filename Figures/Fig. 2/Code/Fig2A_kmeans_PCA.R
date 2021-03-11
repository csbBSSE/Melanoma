#clus_pca function generates a PCA plot marking the pre-defined clusters
#data = dataset containing genes in rows and samples in columns. 
#clus1 = sample list that contains list of samples that fall into the first cluster (based on k-means clustering for top 3000 
#most variable genes for each dataset.)
#title = Title for each plot

clus_pca<- function(data,clus1, title){
  library(ggpubr)
  cut <- is.finite(match(names(data), clus1))   #Assigning cluster number to each sample.
  cut <- replace(cut, cut==T, 1)
  cut <- replace(cut, cut==F, 2)
  
  var_genes=apply(data,1,var)
  data <- data[rev(order(var_genes))[1:3000],]  #top 3000 genes based on variance
  a <- factoextra::fviz_cluster(list(data = t(data), cluster = cut), geom = "point", 
                                pointsize =1.5,  xlab = "PC1", ylab = "PC2", main = title)+
    border()
  return(a)
  
}

ds <- c("GSE7127", "GSE80829", "GSE137391", "CCLE","GSE10916", "GSE4843") 
p<- list()

jpeg("Figures/Fig. 2/2A_Clustering_K_means.jpeg", height = 1200, width = 220)
for ( i in ds){
  
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
  clus1 <- as.character(read.delim(paste0("Datasets/Clusters/", i,".txt"))$x)
  k = which(ds==i)
  p[[k]] <- clus_pca(df,clus1, title = i)
}
cowplot::plot_grid(plotlist = p, ncol = 1)
dev.off()
