#' Function to output the standard error of mean from input vector
#'
#' @param x
#' @keywords data_analysis
#' @export
#' @examples
#' sem()

sem <- function(x) sd(x)/sqrt(length(x))