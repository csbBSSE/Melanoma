# a = dataframe, with genes in column, samples in rows 

bimodality_coefficient <- function(x, na.rm=FALSE) {  #Adapted from modes package
  
  # Remove missing values, if desired
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
  n <- length(x)
  if (n == 0) {
    # The coefficient is not defined for
    # empty data
    return(NaN)  
  } else {
    
    m3 <- psych::skew(x, type=2)
    m4 <- psych::kurtosi(x, type=2)
    
    # Calculate the coefficient based on
    # the above
    bc <- (m3^2 + 1) /
      (m4 + 3 * ((n - 1)^2 / ((n - 2) * (n - 3))))
    
    return(bc)
  }
}

Bimodality <- function(a,list){
  
  require(diptest)
  a <- as.data.frame(t(a))  #Genes in rows, samples in columns
  z <- ncol(a)
  k <- nrow(a)
  n <- data.frame()
  a<- as.matrix(a)
  
  for (i in 1:k) {
    
    x <- a[i,]
    
    #Hartigan's diptest
    m= dip.test(x)
    n[i,1] <- m$statistic
    n[i,2] <- m$p.value                      # n hs values of the statistic, p-value and l is number or non zero values
    
    if (n[i,2] < 0.05){
      
      n[i,3] <- 1
      
    } else {
      
      n[i,3] <- 0
      
    }
    
    #Bimdality coefficient
    n[i,4] <- bimodality_coefficient(x)
    
    
    if (n[i,4]>0.555){
      
      n[i,5] <- 1
      
    } else {
      
      n[i,5] <- 0
      
    }
    
    
    
  }
  rownames(n) <- rownames(a)
  names(n) <- c("Hartigan's dip statistic", "p-value", "Dip test", "BC", "Bimodality Coefficient")
  n <- n[list,]
  ggcorrplot::ggcorrplot(t(n[,c(3,5)]), colors = c("blue", "red", "darkgreen"), show.legend = F)
  
  
  
}

#4C
data <- read.delim("Datasets/RACIPE.txt")
list <- c( "KLF4", "SMAD3","NR3C1", "NR2F1","MITF","TBX3", "NFIC", 
           "TFE3","ETV5",  "FOS",  "JUN","AHR","SMAD4",  "FOXF1", 
           "MAFB",  "STAT5A", "TFAP2A")
jpeg("Figures/Fig. 4/4C_bimodality_coefficient.jpeg")
Bimodality(data,list)
dev.off()

#S4D
data <- as.data.frame(log2(t(read.delim("Datasets/GSE81383.txt"))+1))
list <- c("MITF", "FOS", "ETV5", "SMAD3", "NR2F1","NFIC", "KLF4","JUN", "TFE3", "NR3C1", "TBX3")
data <- data[,list]
jpeg("Figures/Fig. 4/S4D_bimodality_coefficient.jpeg")
Bimodality(data,list)
dev.off()
