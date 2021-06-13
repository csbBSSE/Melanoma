library(monocle)
monocle_sc <- function(df, clusters){
  
  df <- as.data.frame(df)
  clusters <- as.data.frame(clusters)
  rownames(clusters) <- names(df)
  names(clusters) <- c("Phenotype")
  fd <- as.data.frame(rownames(df))
  rownames(fd) <- fd$`rownames(df)`
  names(fd) <- c("gene_short_name")
  pd <- new("AnnotatedDataFrame", data = clusters)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(as.matrix(df), phenoData = pd,featureData = fd, expressionFamily=negbinomial.size())
  
  cds <-  estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds<- detectGenes(cds, min_expr = 0.1)
  ordering_genes <-  c("NFIC", "AHR", "JUN", "SMAD3", "KLF4", "TBX3", "NR3C1","MITF","FOS", "SMAD4","TFAP2A", "NR2F1", "MAFB", "ETV5", "STAT5A",
                       "FOXF1","TFE3") 
  cds <- setOrderingFilter(cds, ordering_genes)
  cds <- reduceDimension(cds, max_components  = 2,
                         method = 'DDRTree', auto_param_selection = F)
  cds <- orderCells(cds)
  cds$Phenotype <- as.factor(cds$Phenotype)
  jpeg("Figures/Fig. 6/6Ci_monocle_phenotype.jpeg")
  print(plot_cell_trajectory(cds, color_by = "Phenotype", show_branch_points = F))
  dev.off()
  jpeg("Figures/Fig. 6/6Cii_monocle_pseudotime.jpeg")
  print(plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F))
  dev.off()
  print(plot_cell_trajectory(cds, color_by = "State"))
  
}

df <- read.delim("Datasets/GSE134432.txt")
clusters <- read.delim("Datasets/Clusters/GSE134432.txt")
monocle_sc(df,clusters)