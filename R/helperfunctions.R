#' Helper functions for stuff...
#'
#' @param various
#' @keywords various
#' @export
#' @examples
nth_highest_val <- function(my_vector,n){ # function to take the n-th highest value of a vector
  my_vec_sorted <- sort(my_vector, decreasing = TRUE)
  return(my_vec_sorted[n])
}

not_in <- function(dataframe_whole,dataframe_part){ # function to extract all parts of dataframe_whole EXCEPT dataframe_part
  val <- subset(dataframe_whole, !(rownames(dataframe_whole) %in% rownames(dataframe_part)))
  return(val) # This needs adjustment for characters...
}

min_max <- function(x,a,b) { # function to normalize data (good in apply of rows)
  if (missing(a))
  {
    a = -1
  }
  if (missing(b))
  {
    b = 1
  }
  mm <- (x - min(x)) / (max(x) - min(x))
  val <- (b-a)*(mm)+a
  return (val)
}

standardize <- function(data_frame){
  data_frame2 <- data_frame
  mean_vals <- apply(data_frame,1,function(x) mean(x,na.rm=T))
  sd_vals <- apply(data_frame,1,function(x) sd(x,na.rm=T))
  for (i in 1:nrow(data_frame)){
    data_frame2[i,] <- lapply(data_frame[i,],function(x) ((x-mean_vals[i])/sd_vals[i]))
  }
  return(data_frame2)
}
