#PCA_hist function generates the PCA plot for all datapoints with density distribution along PC1
#data = dataframe of gene expression, genes in columns
#clus = no. of clusters to be generated
#seed = initialization seed to be set
#pointsize = size of point on PCA plot
#V = Loading scores of genes for principal components

PCA_hist <- function(data,pointsize = 0.5, V = NULL, title){
  
  require(factoextra)
  require(ggpubr)
  
  data[!is.finite(as.matrix(data))] <- 0
  if(is.null(V)){
    
    pca <- FactoMineR::PCA(data, scale = FALSE, graph = F)
    df <- as.data.frame(pca$ind$coord[,1:2])
    
    
  }else{
    
    #Determining coordinates for fixed PC1 loading scores (For control network)
    data <- as.matrix(data)
    df <- data %*% V                                              
    df <- as.data.frame(df)[,1:2]
    
  }
  
  names(df) <- c("PC1", "PC2")
  df <- as.data.frame(df[complete.cases(df),])
  
  sp <- ggscatter(df, x = "PC1", y = "PC2",
                  color = "#0073C299",
                  size = pointsize)+
    border()
  
  sp <- sp + rremove("legend")
  xplot <- ggdensity(df[,1],xlab="PC1", fill = "black")+ ylim(0,0.5)
  
  library(cowplot)
  p <- print(plot_grid(xplot, sp,  ncol = 1, align = "hv", 
                  rel_widths = c(2, 1), rel_heights = c(1, 2)))
  
  title1 <- ggdraw() + draw_label(paste0(title), fontface='bold')
  print(plot_grid(title1, p, ncol=1, rel_heights=c(0.1, 1)) )
  
}

#To determine V (PC1 loading scores)

df <- read.delim("Datasets/KD/Control.txt")
pca <- FactoMineR::PCA(df, graph = F, scale.unit = F)
V <- pca$svd$V

jpeg("Figures/Fig. 7/S6_Control.jpeg", width = 350, height =350)
PCA_hist(df, V = V, title = "Control")
dev.off()

list <- c("NFIC", "AHR", "JUN", "SMAD3", "KLF4", "TBX3", "NR3C1","MITF","FOS", "SMAD4")

for (i in list) {
  
  df <- read.delim(paste0("Datasets/KD/",i," KD.txt"))
  jpeg(paste0("Figures/Fig. 7/S6_",i,"_KD.jpeg"), width = 350, height = 350)
  PCA_hist(df, V = V, title = paste0(i, " KD"))
  dev.off()
  
}


