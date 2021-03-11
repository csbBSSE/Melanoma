#The function corr_circle generates 2 plots - Correlation circle with squared cosine values, and Contribution of each gene to PC1
#data = dataframe for simulated data (genes in columns)
corr_circle <- function(data){
  
  require(FactoMineR)
  require(factoextra)
  require(ggpubr)
  pca <- PCA(data, graph = F)
  jpeg("Figures/Fig. 4/S4A_ii_Corr_Circle.jpeg", width = 450, height = 450)
  print(fviz_pca_var(pca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),            #PCA_correlation circle - S4A_ii
               repel = T, xlab = "PC1", ylab = "PC2", label.size = 0 )+ 
    theme(text = element_text(size = 8),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 8)))
  dev.off()
  
  a <- fviz_contrib(pca, choice = "var", axes = 1, top = 17)                   #Contribution - 4B_iii
  jpeg("Figures/Fig. 4/4B_iii_Cont_PCA.jpeg", width = 600, height = 400)
  print(ggpar(a, font.tickslab =c(16), font.y = c(16)))
  dev.off()
}

data <- read.delim("Datasets/RACIPE.txt")


corr_circle(data)

