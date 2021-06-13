hist <- function(data, genes){
  
  data <- data[, genes]
  data <- as.data.frame(t(data))
  z <- ncol(data)
  k <- nrow(data)
  n <- data.frame()
  #par(mfrow = par)
  for (i in 1:k) {
    
    x <- vector()
    l=1
    for (j in 1:z)  {
      
      zn <- (data[i,j])
      if(zn != 0){
        
        x[l] <- zn    #Add non-zero values to a vector
        l=l+1
        
      }
      
    } 
    print(graphics::hist(x, main = rownames(data)[i], breaks = 25, col = "black", xlab = "levels of expression" , ylab = "Frequency", cex.lab=1.5))
  
  }
  
}  

data <- read.delim("Datasets/RACIPE.txt")


#4C

list <- c("MITF", "FOS", "ETV5", "SMAD3", "NR2F1","NFIC", "KLF4","JUN", "TFE3", "NR3C1", "TBX3") 
jpeg("Figures/Fig. 4/4D_Histogram.jpeg", width = 450, height = 1000)
par(mfrow = c(6,2))
hist(data, list)
dev.off()

#S4A
jpeg("Figures/Fig. 4/S4B_Histogram.jpeg", width = 450, height = 300)
list <- c( "MAFB", "STAT5A", "SMAD4", "FOXF1","TFAP2A","AHR")
par(mfrow = c(2,3))
hist(data, list)
dev.off()

#S4C
data <- as.data.frame(log2(t(read.delim("Datasets/GSE81383.txt"))+1))
list <- c("MITF", "FOS", "ETV5", "SMAD3", "NR2F1","NFIC", "KLF4","JUN", "TFE3", "NR3C1", "TBX3")
jpeg("Figures/Fig. 4/S4C_Histogram_GSE81383.jpeg", width = 500, height = 400)
par(mfrow = c(3,4))
hist(data,list)
dev.off()
