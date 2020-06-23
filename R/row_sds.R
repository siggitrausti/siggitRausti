#' Find row standard deviations
#'
#' @param dataset
#' @keywords data_processing
#' @export
#' @examples
#' row_sds()

row_sds <- function(dataset){
  # A function to find the SDs of rows within a dataset/dataframe/matrix
  outp <- rep(NA,nrow(dataset))
  for (i in 1:nrow(dataset)){
    outp[i] <- sd(dataset[i,],na.rm=T)
  }
  return(outp)
}