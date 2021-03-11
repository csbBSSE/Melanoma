#AIC_BIC function generates a plot for difference for AIC and BIC scores between nth and (n+1)th models for clustering.
#data = dataset containing genes in rows and samples in columns. 
#title = Title for each plot

AIC_BIC <- function(data,title){
  
  var_genes=apply(data,1,var)
  data <- data[rev(order(var_genes))[1:3000],]  #top 3000 genes based on variance
  data <- as.data.frame(t(data))
  data <- scale(data)                           #Standardizing
  
  kmeansAIC = function(fit){      #calculates AIC and BIC for each clustering model
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    return(data.frame(AIC = D + 2*m*k, BIC = D + log(n)*m*k))  
  }
  
  k <-  data.frame()
  for (i in 1:15) {
    
    set.seed(6)
    km <- kmeans(data,i)
    k[i,1:2] <- kmeansAIC(km)  
    
  }

  #Plot difference
  
  m <- data.frame()
  for (i in 1:15){
    
    m[i,1] <-k[i,1] - k[i+1,1]
    m[i,2] <-k[i,2] - k[i+1,2]
    
  }
  m <- m[1:14,]
  rownames(m) <- c(2:15)
  return(plot(seq(2,15),m[,1],xlab="Number of clusters",ylab="Difference",pch=20,cex=2, 
              type=c("l"), main = paste0(title),cex.lab =1.5, ylim = c(-1.5*(10^4), 1.5*(10^4)))+
          lines(seq(2,15),m[,2],pch=20,cex=2, col = 2))
  
  
}

ds <- c("GSE7127", "GSE80829", "GSE137391", "CCLE","GSE10916", "GSE4843")
p <- list() 

jpeg("Figures/Fig. 1/S1B_AIC_BIC", width = 700, height = 500)
par(mfrow = c(2,3)) 
for ( i in ds){
  
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
  k = which(ds==i)

  AIC_BIC(df, title = i)
  
}
dev.off()
