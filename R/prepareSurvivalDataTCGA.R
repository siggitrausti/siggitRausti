#' Preparing data for survival analysis
#'
#' A function to prepare data (genes of interest) for survival analysis.
#' @param dataset,genes,quantile_used,percentage_patients
#' @keywords data_processing
#' @export
#' @examples
#' prepareSurvivalDataTCGA()

prepareSurvivalDataTCGA <- function(dataset,genes,quantile_used,percentage_patients){
  if (colnames(genes) == 'Group1'){
    id_used_3 <- which(rownames(dataset) %in% genes$Group1)
  } else if (colnames(genes) == 'Group2'){
    id_used_3 <- which(rownames(dataset) %in% genes$Group2)
  } else {
    colnames(genes) = 'Group1'
    id_used_3 <- which(rownames(dataset) %in% genes$Group1)
  }
  survival_cl3_data <- dataset[id_used_3,]
  CL3_z <- standardizeTCGA(survival_cl3_data)
  patients_assignment_vector <- rep(0,ncol(CL3_z))
  q_vals <- rep(0,nrow(CL3_z))
  for (j in 1:nrow(CL3_z)){
    q_vals[j] <- quantile(as.numeric(CL3_z[j,3:ncol(CL3_z)]),quantile_used)
  }
  for (i in 3:ncol(CL3_z)){
    q_comp = rep(0,nrow(CL3_z))
    for (j in 1:nrow(CL3_z)){
      q_comp[j] <- CL3_z[j,i] > q_vals[j]
    }
    q_comp <- as.numeric(q_comp)
    if (length(which(q_comp == 1))/nrow(CL3_z) > percentage_patients){
      patients_assignment_vector[i] = 1
    } else {
      patients_assignment_vector[i] = 2
    }
  }
  return(patients_assignment_vector)
}
