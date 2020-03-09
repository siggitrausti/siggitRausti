#' Function to normalize data (good in 'apply' of rows)
#'
#' @param x,a,b
#' @keywords Normalization
#' @export
#' @examples
#' min_max(x,a,b)

min_max <- function(x,a,b) { 
  if (missing(a))
  {
    a = -1
  }
  if (missing(b))
  {
    b = 1
  }
  mm <- (x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T))
  val <- (b-a)*(mm)+a
  return (val)
}