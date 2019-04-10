#' Function to extract the n-th highest value in a vector
#'
#' @param my_vector,n
#' @keywords indexing
#' @export
#' @examples
#' nth_highest_val(my_vector,n)

nth_highest_val <- function(my_vector,n){ # function to take the n-th highest value of a vector
  my_vec_sorted <- sort(my_vector, decreasing = TRUE)
  return(my_vec_sorted[n])
}
