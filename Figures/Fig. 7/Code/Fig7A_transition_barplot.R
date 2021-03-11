#This code generates a barplot for an input of % contribution of each state
#To determine the relative proportions, k-means clustering and heatmaps were used to assign clusters
#The files Control.txt and MITF KD.txt contain the relative proportions of the 4 phenotypes in RACIPE results 


library(dplyr)
library(rstatix)
library(ggpubr)

all <- as.data.frame(read.table("Datasets/Control.txt", header = T))
kd <- as.data.frame(read.table("Datasets/MITF KD.txt", header = T))

#Summarize data
summary <- function(df, net){
  
  summ <- df
  #summ[,4] <- c("U","N","I","M")
  names(summ) <- c("R1", "R2", "R3", "cluster")
  summ <- data.frame(summ$cluster, stack(summ[,1:3]))
  summ[,3] <- net
  names(summ)[3] <- "Network"
  names(summ)[1] <- "cluster"
  return(summ)
}


df <- rbind(summary(all, net = c("Control Network")), summary(kd, net = c("MITF KD")))
df$cluster <- factor(df$cluster, levels = c("M","U","T", "N"))

#Statistical test between the groups
stat.test <- df %>%
  group_by(cluster) %>%
  t_test(values ~ Network)%>%
  add_significance("p")

stat.test <- stat.test %>%
  add_xy_position(x = "cluster", dodge = 0.8)

jpeg("Figures/Fig. 7/7A_transition_barplot.jpeg", width = 700, height = 500)
ggbarplot(
  df, x = "cluster", y = "values", add = "mean_sd", 
  fill= "Network", palette = c("#00AFBB", "#E7B800"),
  position = position_dodge(0.8)) +
  xlab("Phenotype")+
  ylab("Percentage")+
  stat_pvalue_manual(stat.test, x = "cluster",
    label = "p.signif" , tip.length = 0,
    bracket.size = 5, size = 10 )+
  border()+
  theme(text = element_text(size=20))
dev.off()
