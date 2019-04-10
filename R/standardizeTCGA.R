#' Standardizing TCGA data
#'
#' A function to standardize data (x-mean/std)
#' @param data
#' @keywords data_processing
#' @export
#' @examples
#' standardizeTCGA()

standardizeTCGA <- function(data){
  # Remove Na genes:
  data1<- na.omit(data)
  data2 <- data1[,3:ncol(data1)]
  data2 <- replace(data2, data2==0, 0.01)
  data2 <- log2(data2) + 1000
  data3 <- data1[,3:ncol(data1)]
  for (j in 1:nrow(data2)){
    mean_val <- mean(as.numeric(as.character(data2[j,])),na.rm=T)
    std_val = sd(as.numeric(as.character(data2[j,])),na.rm=T)
    for (i in 1:ncol(data2)){
      data3[j,i] = (data2[j,i]-mean_val)/std_val # Possible to increase speed by using apply.. AT least this works.
    }
  }
  data_outp <- cbind(data1[,1:2],data3)
  return(data_outp)
}
