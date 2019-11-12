#' Function to normalize a vector so that 1) No values are below tol_value and 2) total sum of vector is 1 
#'
#' @param iso_vector,tol_value
#' @keywords data_processing
#' @export
#' @examples
#' Raw2MDV()

Raw2MDV <- function(iso_vector,tol_value){
  if (missing(tol_value))
  {
    tol_value <- 0.0001
  }
  for (i in 1:length(iso_vector)){
    if (iso_vector[i] < tol_value){
      iso_vector[i] <- 0
    }
  }
  sum_vec <- sum(iso_vector)
  for (i in 1:length(iso_vector)){
    iso_vector[i] <- iso_vector[i]/sum_vec
  }
  return(iso_vector)
}