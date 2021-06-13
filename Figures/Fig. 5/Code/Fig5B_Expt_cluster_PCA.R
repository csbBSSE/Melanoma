#The function clus_pca generates a PCA plot for a clustered dataset 
#data = sequencing data, genes in rows and samples in columns
clus_pca <- function(data,clus1,genes, title){
  library(ggpubr)
  cut <- is.finite(match(names(data), clus1))
  cut <- replace(cut, cut==T, 1)
  cut <- replace(cut, cut==F, 2)
  
  data <- data[genes,]
  a <- factoextra::fviz_cluster(list(data = t(data), cluster = cut), geom = "point", ellipse = F, xlab = "PC1", ylab = "PC2", main = title)+
    border()
  return(a)
  
}

p <- list()
ds <- c("GSE7127", "GSE80829", "GSE137391", "CCLE","GSE10916", "GSE4843")
list <- c("MITF", "FOS", "ETV5", "SMAD3", "NR2F1","NFIC", "KLF4","JUN", "TFE3", "NR3C1", "TBX3")

for ( i in ds){
  
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1, check.names = F)
  clus1 <- as.character(read.delim(paste0("Datasets/Clusters/", i,".txt"))$x)
  k = which(ds==i)
  p[[k]] <- clus_pca(df,clus1,list, title = i)
  
}

jpeg(filename="Figures/Fig. 5/5B_PCA_2clusters_expt.jpeg", height = 400, width = 600)
cowplot::plot_grid(plotlist = p, nrow = 2)
dev.off()
