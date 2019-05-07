#' Function to standardize dataframe (potentially more global than 'standardizeTCGA)
#'
#' @param data_frame
#' @keywords Normalization
#' @export
#' @examples

standardize <- function(data_frame){
  mean_vals <- apply(data_frame,1,function(x) mean(x,na.rm=T))
  sd_vals <- apply(data_frame,1,function(x) sd(x,na.rm=T))
  if (any(sd_vals == 0)){
    id_remove <- which(sd_vals == 0)
    data_frame <- data_frame[-c(id_remove),]
    mean_vals <- mean_vals[-c(id_remove)]
    sd_vals <- sd_vals[-c(id_remove)]
  }
  data_frame2 <- data_frame
  for (i in 1:nrow(data_frame)){
    data_frame2[i,] <- lapply(data_frame[i,],function(x) ((x-mean_vals[i])/sd_vals[i]))
  }
  return(data_frame2)
}
