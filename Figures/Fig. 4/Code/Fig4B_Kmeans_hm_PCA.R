#The kmeans_hm_pca generates two clusters, and their corresponding heatmap and PCA plot
#k is the number of clusters
#seed is the seed value to initialize the algorithm
#list is the order of genes for the heatmap


kmeans_hm_pca<- function(data, k, seed=123, list, dist =F , save = T){
  
  library(ggpubr)
  set.seed(seed)
  cut <- kmeans(data,k)
  
  if(save == T){
    
    jpeg("Figures/Fig. 4/4B_i_PCA_RACIPE_2.jpeg", width = 560, height = 460)
    print(factoextra::fviz_cluster(cut, geom = "point", data, ellipse =F,  
                             pointsize =1,  xlab = "PC1", ylab = "PC2", labelsize = 28)+
      border())
    dev.off()
    
    jpeg("Figures/Fig. 4/4B_ii_HM_RACIPE_2.jpeg", width = 520, height = 400)
    print(ComplexHeatmap::Heatmap(data, split = cut$cluster, column_order = list, show_row_dend = F , show_row_names = F,
                            heatmap_legend_param = list(title = "Scale")))
    dev.off()
    
    }else {
      print(factoextra::fviz_cluster(cut, geom = "point", data, ellipse =F,  
                                   pointsize =1,  xlab = "PC1", ylab = "PC2", labelsize = 28)+
              border())
      print(ComplexHeatmap::Heatmap(data, split = cut$cluster, column_order = list, show_row_dend = F , show_row_names = F,
                                  heatmap_legend_param = list(title = "Scale"))) 
    
  }
  
 
 
}


data <- read.delim("Datasets/RACIPE.txt")
list <- c("MITF", "FOS", "SMAD4", "STAT5A", "ETV5",  
          "SMAD3", "NR2F1", "NFIC",  "KLF4", "JUN","TFE3","NR3C1",
          "MAFB",  "TBX3","FOXF1","AHR", "TFAP2A")


kmeans_hm_pca(data,k=2, seed = 44, list =list)
