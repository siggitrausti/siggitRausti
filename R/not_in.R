#' Function to extract all parts of dataframe_whole EXCEPT dataframe_part
#'
#' @param dataframe_whole,dataframe_part
#' @keywords intersect
#' @export
#' @examples
not_in <- function(dataframe_whole,dataframe_part){ 
  val <- subset(dataframe_whole, !(rownames(dataframe_whole) %in% rownames(dataframe_part)))
  return(val) # This needs adjustment for characters...
}