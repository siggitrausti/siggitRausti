#' Installation check of packages
#'
#' A function for providing the necessary packages for the functions it uses...
#' @param package_names
#' @keywords packages
#' @export
#' @examples
#' load_or_install('ggplot2')

load_or_install<-function(package_names){
  for(package_name in package_names)
  {
    if(!is_installed(package_name))
    {
      install.packages(package_name)
    }
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }
}