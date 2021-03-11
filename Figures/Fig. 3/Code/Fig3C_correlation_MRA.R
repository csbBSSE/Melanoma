#correlation function
#df = dataframe  for sequencing data, genes in rows, samples in columns 
#list = list of genes to be plotted
#method = "pearson", "spearman"
#save = TRU/FALSE, to save image 
#z = Renames and replaces last element (used for ZEB1, since it is identified under different annotations)

correlation <- function(df, list, method, save = FALSE, title, z = F){
  require(ggplot2)
  require(ggcorrplot)
  
  df<- df[list,]
  
  if(z==T){
    
    if(rownames(df)[nrow(df)]==list[length(list)]){ #Converts the last element to second last, for ZEB1
      
      df <- df[-c((nrow(df)-1)),]
      rownames(df)[nrow(df)] <- list[length(list)-1]  #Renaming "LOC100996668 /// ZEB1" to "ZEB1"
      
    }
  }
  
  df<- df[complete.cases(df),]  # Removing NAs, in case they exist
  z <- ncol(df)
  k<- nrow(df)
  
  corr <- cor(t(df) , method = method)  #Correlation matrix
  corr.p <- cor_pmat(t(df), method = method ) #p-values
  
  if(save == T){
    
    jpeg(paste0(title,"_corr_",method,".jpeg"))
    ggcorrplot(corr,
               hc.order = FALSE,
               type = "lower",
               outline.color = "white", p.mat = as.matrix(corr.p), sig.level = 0.05, insig = "pch",  #
               show.diag = TRUE, show.legend = FALSE)
    dev.off()
    
  }else{
    
    return(ggcorrplot(corr,
                     hc.order = FALSE,
                     type = "lower",
                     outline.color = "white", p.mat = as.matrix(corr.p), sig.level = 0.05, insig = "pch",
                     show.diag = TRUE, show.legend = FALSE, title = paste0(title," (n=",ncol(df),")")))
    
  }
  
  
}


ds <- c("GSE4843","GSE137391")    
list <- as.character(read.delim("Datasets/MRA_genes.txt", header = T)[,1])

p <- list()

for ( i in ds){
  
  k = which(ds==i)
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
  p[[k]] <- correlation(df, method = "spearman", list = list, title =i, z=F) # use method = "spearman" or "pearson"
  
}

jpeg("Figures/Fig. 3/3C_Correlation_MRA.jpeg",width = 1050, height = 680 )
cowplot::plot_grid(plotlist = p, nrow = 1)
dev.off()

ds <- c("GSE7127","GSE80829", "CCLE", "GSE10916")   
p <- list()

for ( i in ds){
  
  k = which(ds==i)
  df <- read.delim(paste0("Datasets/", i,".txt"), row.names = 1)
  p[[k]] <- correlation(df, method = "spearman", list = list, title =i, z=F) # use method = "spearman" or "pearson"
  
}

jpeg("Figures/Fig. 3/S3B_Correlation_MRA.jpeg",width = 1050, height = 680 )
cowplot::plot_grid(plotlist = p, nrow = 2)
dev.off()

