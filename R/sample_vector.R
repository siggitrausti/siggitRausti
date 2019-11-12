#' A function for sampling n rows from dataframe. 
#'
#' @param vector,n,x,replicate
#' @keywords data_processing
#' @export
#' @examples
#' sample_vector()

sample_vector <- function(vector,n,x,replicate=T){
  if (missing(replicate))
  {
    replicate <- FALSE
  }
 # vector is the vector of choice, n is the number of samples and x the number of times n should be sampled from vector
  # Creates a n*x dataframe of sampled data
  length_vector <- length(vector)
  outp <- c()
  for (i in 1:x){
    if (replicate == FALSE){
      outp <- cbind(outp,vector[sample(1:length_vector, n,replace = FALSE)])
    } else {
      outp <- cbind(outp,vector[sample(1:length_vector, n,replace = TRUE)])
    }
  }
  return(outp)
}