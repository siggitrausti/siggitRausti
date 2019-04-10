#' Differential gene analysis
#'
#' A function for extracting DE genes from TCGA data. This is assuming the data is normalized and uses
#' students t-test to evaluate, and BH-FDR correction. FC refers to log-fold change...
#' group1 and group2 refer to different cluster of samples from same dataframe. 'return_val' tells you
#' whether to take the significantly higher genes from group1 or group2 (1 and 2, respectively).
#' @param group1,group2,fc,pval,return_val
#' @keywords DEgenes
#' @export
#' @examples
#' extractDEgenesTCGA()

extractDEgenesTCGA <- function(group1,group2,fc,pval,return_val){

  if (missing(return_val))
  {
    return_val <- 1
  }
  if (dim(group1)[1] == 0){
    stop("No genes involved with input group!")
  }
  novar1<- which(apply(group1, 2, var) == 0)
  novar2<- which(apply(group2, 2, var) == 0)
  novar_fin = union(novar1,novar2)
  group1[,novar_fin] = NULL # to def take out rows with zero variance...
  group2[,novar_fin] = NULL # -||-

  group1_2 <- replace(group1, group1==0, 0.01)
  group2_2 <- replace(group2, group2==0, 0.01)
  log_dat1 <- log2(group1_2) + 1000
  log_dat2 <- log2(group2_2) + 1000
  #print(t.test(log_dat1[,2],log_dat2[,2],alternative="two.sided",var.equal=F))
  p_values = rep(0,ncol(group1_2))
  for (i in 1:ncol(group1_2)){
    ttest_res <- t.test(log_dat1[,i], log_dat2[,i], alternative = "two.sided", var.equal = FALSE) # Should the variance be equal?
    p_values[i] <- ttest_res$p.value
  }
  # Correct p-values for multiple comparisons (Benjamini-Hochberg...)
  adjusted_pvals <- round(p.adjust(p_values,'BH',length(p_values)),3)

  # Find the genes that are at least twofold higher in cluster 1 and 2 and are significant:
  sign_in_cluster1 <- c()
  sign_in_cluster2 <- c()
  for (i in 1:ncol(group1_2)){
    log_dat_diff <- mean(log_dat1[,i])-mean(log_dat2[,i])
    if (log_dat_diff >= fc && adjusted_pvals[i] < pval){
      sign_in_cluster1 <- c(sign_in_cluster1,colnames(group1_2)[i])
    } else if (log_dat_diff <= -fc && adjusted_pvals[i] < pval){
      sign_in_cluster2 <- c(sign_in_cluster2,colnames(group2_2)[i])
    }
  }
  cluster1 <- data.frame(sign_in_cluster1)
  if (dim(cluster1)[1] == 0){
    print('No significant genes!')
    return(data.frame(NULL))
  } else {
    colnames(cluster1) <- 'Group1'
    cluster2 <- data.frame(sign_in_cluster2)
    colnames(cluster2) <- 'Group2'
    if (return_val == 1){
      return(cluster1)
    } else if (return_val == 2){
      return(cluster2)
    }
  }
}
