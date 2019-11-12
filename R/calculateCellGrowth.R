#' Get growth rate from cell numbers and time
#'
#' @param cell_no,timepoints
#' @keywords Cells
#' @export
#' @examples
#' calculateCellGrowth()

calculateCellGrowth <- function(cell_no,timepoints){
  # This function is to calculate growth rate from cell numbers and hours
  if(length(cell_no) != length(timepoints)){
    print('No of cells and timepoints dont match, try again stupid...')
    slope <- c()
  } else {
    GROWTH <- timepoints # This is for the beauty of the output solely....
    log_vals <- log(cell_no)
    slope <- coef(lm(log_vals ~ GROWTH))[2]
  }
  return(slope)
}