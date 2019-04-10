#' Hierarchical clustering of dataframe
#'
#' A function for Hierarchical clustering (Spearman or Pearson) of data from cBioPortal (raw data).
#' Useful for looking at specific patients or all of them.
#' @param dataframe,gene_id,cor_method
#' @keywords clustering
#' @export
#' @examples
#' Hierarchical_clustering()

Hierarchical_clustering <- function(dataframe,gene_id,cor_method){
  load_or_install('cluster')
  data_zscored2 <- na.omit(dataframe)
  if(gene_id == 'Hugo'){
    data_zscored2$Entrez_Gene_Id <- NULL
  }
  else if(gene_id == 'Entrez'){
    data_zscored2$Entrez_Gene_Id <- NULL
  }
  if(cor_method == 'Spearman'){
    data_zscored_3.cor <- cor(data_zscored2[,-1],method="spearman") # Creates the correlation matrix
  }
  else if(cor_method == 'Pearson'){
    data_zscored_3.cor <- cor(data_zscored2[,-1],method="pearson") # Creates the correlation matrix
  }
  d2.cor <- as.dist(1 - data_zscored_3.cor)
  data_zscored_3.hclust <- hclust(d2.cor, method = "ward.D2")
  plot(data_zscored_3.hclust, cex = 0.6,main=paste('Hierarchical Clustering (',cor_method,')'),xlab="",sub="")
  return(data_zscored_3.hclust)
}
