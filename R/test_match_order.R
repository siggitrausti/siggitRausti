#' Function to check whether two vectors are in the same order or not. 
#'
#' @param x,y
#' @keywords data_processing
#' @export
#' @examples
#' test_match_order()

test_match_order <- function(x,y) {
  
  if (all(x==y)) print('Perfect match in same order')
  
  if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')
  
  if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}