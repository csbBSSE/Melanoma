#The function cluster_Expression generates heatmaps for a list of gene, across 
#defined clusters. 
#data = dataframe with samples in columns, genes in rows
#clus1 = text containing names of samples in cluster 1
#genes = list of genes to consider for PCA
#title = title of plot


cluster_Expression <- function(data, clus1, genes,title){
  
  require(grid)
  data <- data[genes,]
  data <- scale(t(data))                                   #Z-normalize data
  cut <- is.finite(match(rownames(data),clus1))         #Assign tags to each row based on cluster
  cut <- as.integer(cut)
  cut <- replace(cut, cut==1,2)
  cut <- replace(cut, cut==0,1)
  print(ComplexHeatmap::Heatmap(as.matrix(data), split = cut, show_row_names=FALSE, column_order = genes, 
                                heatmap_legend_param = list(title = "Scale"), column_title_side = "top",
                                column_title = paste0(title, " (n = ", nrow(data),")"), 
                                column_title_gp = gpar(fontsize = 24),
                                column_names_gp = gpar(fontsize = 24), show_row_dend = F))
  
}


#5B
ds <- c("GSE4843","GSE137391")
list <- c("MITF", "FOS", "ETV5", "SMAD3", "NR2F1","NFIC", "KLF4","JUN", "TFE3", "NR3C1", "TBX3")

for ( i in ds){
  
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
  clus1 <- as.character(read.delim(paste0("Datasets/Clusters/", i,".txt"))$x)
  k = which(ds==i)
  jpeg(filename=paste0("Figures/Fig. 5/5A_cluster_expr_",i,".jpeg"), height = 500, width = 700)
  cluster_Expression(df,clus1,list, i)
  dev.off()
  
}



#S4D
ds <- c("GSE7127", "GSE80829")
list <- c("MITF", "FOS", "ETV5", "SMAD3", "NR2F1","NFIC", "KLF4","JUN", "TFE3", "NR3C1", "TBX3")
p <- list()

for ( i in ds){
  
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
  clus1 <- as.character(read.delim(paste0("Datasets/Clusters/", i,".txt"))$x)
  k = which(ds==i)
  jpeg(filename=paste0("Figures/Fig. 5/S4C_cluster_expr_",i,".jpeg"), height = 500, width = 700)
  cluster_Expression(df,clus1,list, i)
  dev.off()
  
  
}

