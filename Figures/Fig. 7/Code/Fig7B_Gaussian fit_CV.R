#The following code fits gaussian curves to distribution of points along PC1. 
#Based on the fits, the Coefficient of variance for the two main peaks is estimated and compared

library(dplyr)
library(rstatix)
library(ggpubr)

#To determine PCA loadings (V)
df <- read.delim("Datasets/KD/Control.txt")
pca <- FactoMineR::PCA(df, graph = F, scale.unit = F)
V <- pca$svd$V

#Fitting Gaussian and estimatinf CV
ds <- c("Control", "MITF KD")
k <- c(3,2)
CV <- data.frame()
rep <- c("", "_rep_2","_rep_3")

for (i in ds) {
  
  library(mixtools)
  data <- read.delim(paste0("Datasets/KD/", i,".txt"))
  data[!is.finite(as.matrix(data))] <- 0
  df <- as.matrix(data)%*%V
  df <- as.data.frame(df)[,1]
  set.seed(20)
  my_mix = normalmixEM(df, k=k[which(ds==i)])
  
  jpeg(paste0("Figures/Fig. 7/7B_dist_",i,".jpeg")) #Fit is saved only for 1 replicate
  plot(my_mix,which=2, cex.axis = 1.2, cex.lab = 1.2)
  dev.off()
  
  for (r in rep){
    
    path <- paste0("Datasets/KD/", i,r,".txt")
    r_1 <- which(rep == r)
    data <- read.delim(path)
    data[!is.finite(as.matrix(data))] <- 0
    df <- as.matrix(data)%*%V
    df <- as.data.frame(df)[,1]
    set.seed(20)
    my_mix = normalmixEM(df, k=k[which(ds==i)])    #For fitting, k=3 and k=2 is taken for optimal fitting of Control and MITF KD, respectively
    plot(my_mix,which=2, cex.axis = 1.2, cex.lab = 1.2)                       #To plot fit of all datasets
    
    #Generation of dataframe containing CV values
    
    if (i == "Control"){    
      
      CV[r_1,1:2] <- my_mix$sigma[c(3,1)]/my_mix$mu[c(3,1)] 
      
    }else{
      
      CV[3+r_1,1:2] <- my_mix$sigma/my_mix$mu
      
    }
    
  }
 
}


CV[1:3,3] <- c("Control")
CV[4:6,3] <- c("MITF KD")
names(CV) <- c("Peak 1", "Peak 2", "Network")
a <- reshape2::melt(CV)
a$value <- abs(a$value)
a$variable <- factor(a$variable, levels = c("Peak 1","Peak 2"))
names(a)[2] <- c("Peak")

#Statistical testing
stat.test <- a%>%
  group_by(Peak) %>%
  t_test(value ~ Network)%>%
  add_significance("p")

stat.test <- stat.test %>%
  add_xy_position(x = "Peak", dodge = 0.8)


p <- ggbarplot(data = a, x= "Peak", y="value",fill= "Network",add = "mean_sd",
               palette = c("#00AFBB", "#E7B800"),position = position_dodge(0.9))+
  stat_pvalue_manual(stat.test, x = "Peak",
                   label = "p.signif" , size = 10 )+
  border()

jpeg("Figures/Fig. 7/7B_Barplot.jpeg", width = 560, height = 475)
p +rremove("legend.title") + rremove("x.title") + ylab("Absolute Coefficient of Variance")
dev.off()



