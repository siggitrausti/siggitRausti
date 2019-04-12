#' Installation check of packages
#'
#' This function allows you to check whether a certain package is installed or not
#' @param mypkg
#' @keywords packages
#' @export
#' @examples
#' is_installed('ggplot2')

is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])