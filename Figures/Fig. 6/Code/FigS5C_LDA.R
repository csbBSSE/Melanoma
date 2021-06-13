#This code generates a plot of absolute value of contribution of each gene
#to LD1
#df = dataframe containing samples to be compared, samples are rows, genes are columns
#cut = List of cluster symbols for each sample in df

library(dplyr)
library(ComplexHeatmap)


lda_val <- function(df,cut){
  
  lda <- MASS::lda(df, grouping = as.factor(cut))
  d <- as.data.frame(lda$scaling)
  print(barplot(abs(d$LD1), col = "coral", ylab = "Absolute value of loading score",las=2, names.arg = rownames(d),
                ylim=range(pretty(c(0, abs(d$LD1)))), cex.lab = 1.5, cex.names = 1.2))
  return(d)
}

df <- read.delim("Datasets/RACIPE.txt")
list <- c( "KLF4", "SMAD3","NR3C1", "NR2F1","MITF","TBX3", "NFIC", 
           "TFE3","ETV5",  "FOS",  "JUN","AHR","SMAD4",  "FOXF1", 
           "MAFB",  "STAT5A", "TFAP2A")

set.seed(44)
cut <- kmeans(df, 4)$cluster   #Separation of clusters
data1 <- df[cut==1|cut==2,]    #Clubbing proliferative subclusters
data2 <- df[cut==4|cut==3,]    #Clubbing invasive subclusters

cut1 <- cut[cut==1|cut==2]
cut2 <- cut[cut==4|cut==3]

jpeg("Figures/Fig. 6/S5C_LDA_Pro.jpeg", width = 830, height = 480)
lda_val(data1[,list],cut1)
box()
abline(h=0.6)
dev.off()

jpeg("Figures/Fig. 6/S5C_LDA_Inv.jpeg", width = 830, height = 480)
lda_val(data2[,list],cut2)
box()
abline(h=0.6)
dev.off()