silhouette <- function(data, range =c(2:10)){
  
  d = dist(data)^2
  a <- data.frame()
  
  for (k in range){
    for (j in 1:3){
      km.res <- kmeans(data, centers = k)
      ss <- cluster::silhouette(km.res$cluster, d)
      a[j,k-1] <- mean(ss[,3])
      
    }
    
  }
  names(a) <- range
  print(boxplot(a, main = " Silhouette width plot", xlab = "K- value", ylab = "Silhouette width") )
  
}

data <- read.delim("Datasets/RACIPE.txt")
jpeg("Figures/Fig. 4/S4A_i_Silhouette.jpeg")
silhouette(data)
dev.off()