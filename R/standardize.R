#' Function to standardize dataframe (potentially more global than 'standardizeTCGA)
#'
#' @param data_frame
#' @keywords Normalization
#' @export
#' @examples

standardize <- function(data_frame){
  data_frame2 <- data_frame
  mean_vals <- apply(data_frame,1,function(x) mean(x,na.rm=T))
  sd_vals <- apply(data_frame,1,function(x) sd(x,na.rm=T))
  if (any(sd_vals == 0)){
    data_frame[which(sd_vals == 0),] = NULL
  }
  for (i in 1:nrow(data_frame)){
    data_frame2[i,] <- lapply(data_frame[i,],function(x) ((x-mean_vals[i])/sd_vals[i]))
  }
  return(data_frame2)
}
