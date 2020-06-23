#' A function to remove rows with zero variance (where sd is zero)
#'
#' @param dataset,threshold,rows_or_cols
#' @keywords data_processing
#' @export
#' @examples
#' remove_low_variance()

remove_low_variance <- function(dataset,threshold,rows_or_cols){
  if (missing(threshold))
  {
    threshold <- 0.1
  }
  if (missing(rows_or_cols))
  {
    rows_or_cols <- 'Rows'
    print('Missing rows_or_cols variable. Assuming its rows...')
  }
  if (rows_or_cols == 'Rows'){
    coef.var <- rep(0,nrow(dataset))
    for (i in 1:nrow(dataset)){
      coef.var[i] = sd(as.numeric(dataset[i,]),na.rm=T)/mean(as.numeric(dataset[i,]),na.rm=T)
    }
    dataset2 <- data.frame(dataset[which(coef.var>threshold),])
    return(dataset2)
  } else if (rows_or_cols == 'Cols'){
    coef.var <- rep(0,ncol(dataset))
    for (i in 1:ncol(dataset)){
      coef.var[i] = sd(as.numeric(dataset[,i]),na.rm=T)/mean(as.numeric(dataset[,i]),na.rm=T)
    }
    dataset2 <- data.frame(dataset[,which(coef.var>threshold)])
    return(dataset2)
  } else {
    return(NULL)
  }
  
}