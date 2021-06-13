library(ggcorrplot)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggplot2)

ds <- c("GSE10916","GSE4843","CCLE", "GSE80829","GSE7127",  "GSE137391" ) 
ds_MRA <- c("GSE4843", "GSE10916", "CCLE")

for (k in ds_MRA) {
  list <- read.delim(paste0("Datasets/MRA genes/",k,".txt"))[,1]
  df <- read.delim(paste0("Datasets/",k,".txt"), row.names = 1)
  df<- df[list,]
  df<- df[complete.cases(df),]  # Removing NAs, in case they exist
  method <- "spearman"
  corr <- as.data.frame(cor(t(df) , method = method))  #Correlation matrix
  n <- length(list)
  prol <- names(corr)[corr$MITF >0] #Set order for genes -use NFATC4 for CCLE
  inv <- names(corr)[!names(corr)%in% prol]
  team <- vector()
  ds_other <- ds[ds!=k]
  for (i in ds_other) {
    df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
    df <- df[list,]
    df <- df[complete.cases(df),]
    method <- "spearman"
    corr <- as.data.frame(cor(t(df) , method = method))  #Correlation matrix
    corr[is.na(corr)] <- 0
    inv_sub <- inv[inv %in% names(corr)]
    prol_sub <-prol[prol %in% names(corr)]
    corr <- as.matrix(corr)
    team[which(ds_other==i)] <- (sum(corr[prol_sub,prol_sub]) + sum(corr[inv_sub,inv_sub]) - sum(corr[inv_sub,prol_sub]))/n^2
    
  }
  names(team) <- ds_other
  write.csv(as.data.frame(team), paste0("Figures/Fig. 3/Team strenghth_",k,".csv"))
}

df <- data.frame(c(1:5))
for (i in ds_MRA) {
  df <- cbind(df,read.csv(paste0("Figures/Fig. 3/Team strenghth_",i,".csv"))[,2])
}
df <- df[,-1]
names(df) <- ds_MRA
df <- reshape2::melt(df)
names(df) <- c("Dataset", "TS")

stat.test <- df %>%
  t_test(TS ~ Dataset)
stat.test <- stat.test %>%
  add_xy_position(x = "Dataset", dodge = 0.8)


ggbarplot(
  df, x = "Dataset" , y ="TS", add = "mean_sd", 
  fill= "Dataset",title = "Team Strength",
  position = position_dodge(0.8)) +
  xlab("Dataset")+
  ylab("Team strength/(Number of genes)^2")+
  stat_pvalue_manual(stat.test,
                     label = "p" , tip.length = 0 )
ggsave(paste0("Figures/Fig. 3/S3C_Difference_Team strength_datasets.jpeg"))
