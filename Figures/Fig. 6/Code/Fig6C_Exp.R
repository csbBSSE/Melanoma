#This code needs manually reading the datasets GSE134432 and assigning cluster phenotypes.
#This code generates the figure or GSE4843
#For GSE134432, clusters as derived in Wouters et al.,2020 and data was CPM normalized 
markers <- c("AXL", "NGFR", "MLANA", "ETV4")                       #Set of genes on the basis of which clustering is done
list <- c("MITF", "TBX3", "NR2F1", "KLF4", "FOXF1")                #Genes to check expression levels

#Clustering based on marker genes - Only used for GSE4843
clus_marker <- function(data, markers, seed =123, plot =T){
  
  set.seed(seed)
  cut <- kmeans(data[,markers],4)
  
  if(plot ==T){
    
    print(ComplexHeatmap::Heatmap(data[,markers], split = cut$cluster, column_order = markers,
                                  show_row_dend = F , heatmap_legend_param = list(title = "Scale"),
                                  show_row_names = F))
    
  }
  
  return(cut)
}

#Expression function
#cut represent the cluster symbol assigned to each phenotype
expr <- function(data, list, cut, title){                        
  
  require(ggpubr)
  data <- data[,list]
  exp <- as.data.frame(cbind(data,cut))
  exp <- data.frame(exp$cut, stack(exp[,list]))
  exp$exp.cut <- factor(exp$exp.cut, levels = c("U","N","M","T"))
  my_comparisons <- list(  c("U","N"),  c("T","M"))
  p <- ggboxplot(exp, x = "exp.cut", y = "values",title = title,
                 fill = "exp.cut", palette = "jco", 
                 line.color = "gray", line.size = 0.4, 
                 facet.by = "ind", short.panel.labs = FALSE, nrow =1, 
                 panel.labs = list(ind = c("MITF", "TBX3", "NR2F1", "KLF4", "FOXF1")))+
    ylab("Z-score")+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", label.y = 7.5,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                          symbols = c("****", "***", "**", "*", "ns")))+
    ylim(c(-2,10))+
    theme(axis.title = element_text(size = 16))
  
  p <- ggpar(p,legend.title = "Phenotype", legend = "right")    
  
  return(p +  rremove("x.text") + theme(strip.text.x = element_text(size = 16)))
}


#PCA
data <- read.delim(file.choose(), row.names = 1)  #Genes in rows and samples in columns
data <- as.data.frame(scale(t(data)))             #Genes in columns, samples in rows
cut <- clus_marker(data,markers,seed = 8)$cluster #For GSE4843, for GSE134432 set pre-assigned cluster labels

cut[cut==1] <- c("M")
cut[cut==2] <- c("T")
cut[cut==3] <- c("N")
cut[cut==4] <- c("U")

jpeg("Figures/Fig. 6/6C_avg_exp_4843.jpeg", width = 950, height = 300)
expr(data, list, cut, "GSE134432") 
dev.off()

