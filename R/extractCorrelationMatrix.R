#' A function to identify genes in a gene expression matrix. The dataset needs to have genes in columns and patients/cells in rows
#'
#' @param dataset
#' @keywords data_processing
#' @export
#' @examples
#' extractCorrelationMatrix()

extractCorrelationMatrix <- function(dataset){
  require('WGCNA')
  source('./removeLowVariance.R')
  dataset <- removeLowVariance(dataset,0.0005,'Cols') # 0.01 here by default...
  cor_mat <- data.frame(WGCNA::cor(dataset,use='pairwise.complete.obs',method='pearson'))
  return((cor_mat)) # NOTE: ADD SO MULTIPLE GENES CAN BE TESTED AT ONCE, AND GENERATE A LIST WITH 
    # GENES CORRELATED WITH EACH GENE OF THOSE MULTPLE GENES. HAVE IT IN HERE SO THAT THE CORRELATION MATRIX (THE MOST TIME-CONSUMING 
    # PROCESS) IS NOT PERFORMED MORE THAN ONCE. IT IS PERHAPS DUMB TO DO THIS, MAYBE JUST EXPORT THE WHOLE CORRELATION MATRIX, 
    # AND INCLUDE AN OPTION ON WHETHER OR NOT SPECIFIC GENES ARE WANTED...
}