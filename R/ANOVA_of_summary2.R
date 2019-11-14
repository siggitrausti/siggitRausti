#' Compare 3 or more groups using ANOVA and only the means + stds
#'
#' @param groups,means,stds,var.equal
#' @keywords data_processing
#' @export
#' @examples
#' ANOVA_of_summary2()

ANOVA_of_summary2 <- function(groups,means,stds,var.equal=T){
  # This function creates three datapoints from mean and std, and performs either ANOVA of Kruskal-Wallis
  require(reshape2)
  if (missing(var.equal)){
    var.equal = T
  }
  data_points <- data.frame(matrix(nrow=length(means),ncol=4))
  data_points[,1] <- groups
  colnames(data_points)[1] <- 'Groups'
  for (i in 1:length(means)){
    data_points[i,2:ncol(data_points)] <- c(means[i]+stds[i],means[i],means[i]-stds[i])
  }
  
  mdata <- melt(data_points, id=c("Groups"))
  colnames(mdata)[2] <- 'Replicates'
  if(var.equal == T){
    res.aov <- aov(value ~ Groups, data = mdata)
    if(summary(res.aov)[[1]][["Pr(>F)"]][1] >= 0.05){
      print('No significance observed in ANOVA. Tukey test not appropriate...')
    } else {
      outp <- TukeyHSD(res.aov)
    }
  } else {
    print('Dont use ANOVA for non-parametric data! Use Kruskal-Wallis...')
  }
  return(outp)
  }