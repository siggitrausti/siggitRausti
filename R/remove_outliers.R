#' Remove outliers from vector above 1.25 x IQR
#'
#' @param x,na.rm
#' @keywords data_processing
#' @export
#' @examples
#' remove_outliers()

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25,.5,.75), na.rm = na.rm, ...)
  median_val <- 
  H <- 1.25 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- qnt[2]
  y[x > (qnt[3] + H)] <- qnt[2]
  return(y)
}
