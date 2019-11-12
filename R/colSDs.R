#' Find column standard deviations
#'
#' @param dataset
#' @keywords data_processing
#' @export
#' @examples
#' colSDs()

colSDs <- function(dataset){
  # A function to find the SDs of rows within a dataset/dataframe/matrix
  outp <- rep(NA,ncol(dataset))
  for (i in 1:ncol(dataset)){
    outp[i] <- sd(dataset[,i],na.rm=T)
  }
  return(outp)
}