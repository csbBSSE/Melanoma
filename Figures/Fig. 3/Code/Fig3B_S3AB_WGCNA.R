library(WGCNA)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ComplexHeatmap)

data <- read.delim("Datasets/GSE4843.txt", row.names = 1)
data <- as.data.frame(t(data))

#sampling genes
gsg = goodSamplesGenes(data, verbose = 3);
data = data[gsg$goodSamples, gsg$goodGenes]


#soft threshold - identified as 4
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(data, powerVector = powers, verbose = 3)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.90,col="red")         #cutoff
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#Module processing
disTOM <-(1-TOMsimilarityFromExpr(data, power = 4, corType = "pearson", TOMType = "unsigned"))   #TOM dissimilarity matrix
tree <- hclust(as.dist(disTOM), method= "average")
col <- cutreeDynamic(dendro = tree,distM =disTOM, cutHeight = 0.995,   #Adaptive pruning of dendrogram
                     deepSplit = 2, pamRespectsDendro = FALSE,
                     minClusterSize =100);
col <- labels2colors(col) 

sizeGrWindow(8,6)
plotDendroAndColors(tree, col, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Melanoma WGCNA")


#Module eigengenes and merging dissimilar modules
MEList = moduleEigengenes(data, colors = col) 
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")      
merge = mergeCloseModules(data, col, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)

jpeg("Figures/Fig. 3/S3A_i_Dendrogram_WGCNA")
plotDendroAndColors(tree,  mergedColors,
                    c("Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()




#Assigning clusters to module eigengenes
clus1 <- as.character(read.delim("Datasets/Clusters/GSE4843.txt")$x)
cut <- is.finite(match(rownames(data), clus1))
cut <- replace(cut, cut==T, 1)
cut <- replace(cut, cut==F, 2)

ME1 <- mergedMEs[cut==1,]
ME2 <- mergedMEs[cut==2,]
ME1 <- cbind(ME1,c("coral2"))
names(ME1)[(ncol(ME1))] <- c("cut")
ME2 <- cbind(ME2,c("cyan3"))
names(ME2)[(ncol(ME2))] <- c("cut")

MEfin <- rbind(ME1,ME2)
n <- ncol(MEfin)

#Selecting relevant modules - Comparing expression between proliferative and invasive samples
names(MEfin) <- gsub("ME","",names(MEfin))
df <- data.frame(MEfin$cut, stack(MEfin[,1:29]))
names(df)[1] <- c("cluster")

stat.test <- df %>%
  group_by(ind) %>%
  t_test(values ~ cluster)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj", cutpoints = c(0, 1e-04, 0.001, 0.01, 1),
                   symbols = c( "***", "**", "*", "ns"))

stat.test <- stat.test %>%
  add_xy_position(x = "ind", dodge = 0.8)

jpeg("Figures/Fig. 3/S3A_ii_Eigengene_comparison.jpeg", width = 1000, height =500)
ggboxplot(
  df, x = "ind", y = "values",  
  fill= "cluster", palette = c("#ee6a50", "#00cdcd"))+
  xlab("Module")+
  ylab("Eigengene")+
  
  stat_pvalue_manual(stat.test,
                     label = "p.adj.signif" , tip.length = 0)+
  border()+
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "right")+
  scale_fill_discrete(name ="Phenotype", labels = c("Proliferative", "Invasive"))
dev.off()


jpeg("Figures/Fig. 3/S3A_iii_Eigengene_scatter.jpeg", width = 500, height =250)
plot(MEfin$yellow,MEfin$salmon, xlab = "Yellow eigengene", ylab = "Salmon eigengene", col= as.character(MEfin$cut), pch = 16,
     ylim = c(-0.2,0.4))
dev.off()

#Heatmap and Module eigengene value
#Salmon module
data_hm <- scale(data[,mergedColors=="salmon"])
data_hm <- data_hm[rownames(MEfin),]
jpeg("Figures/Fig. 3/3B_ii_Salmon_Eigengene.jpeg", width = 500, height =250)
barplot(MEfin[,"salmon"], col = as.character(MEfin$cut))
dev.off()

jpeg("Figures/Fig. 3/3B_ii_Salmon_Heatmap.jpeg", width = 500, height =250)
Heatmap(t(data_hm), column_split =MEfin$cut, show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, show_column_names = FALSE,
        show_parent_dend_line = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, column_title = NULL, 
        heatmap_legend_param = list(title= c("Scale")))
dev.off()

#Yellow module
data_hm <- scale(data[,mergedColors=="yellow"])
data_hm <- data_hm[rownames(MEfin),]
jpeg("Figures/Fig. 3/3B_i_Yellow_Eigengene.jpeg", width = 500, height =250)
barplot(MEfin[,"yellow"], col = as.character(MEfin$cut))
dev.off()

jpeg("Figures/Fig. 3/3B_i_Yellow_Heatmap.jpeg", width = 500, height =250)
Heatmap(t(data_hm), column_split =MEfin$cut, show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, show_column_names = FALSE,
        show_parent_dend_line = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,  column_title = NULL, 
        heatmap_legend_param = list(title= c("Scale")))
dev.off()




