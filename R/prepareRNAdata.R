#' Preparing RNA data for PCA plotting
#'
#' Function for preparing the RNA data for metPCA acc. to protocol.
#' rna_data' is the raw RNAseq data from cBioPortal. Clinical data is the corresponding information on patients.
#' groups is how many groups the patients should be split into based on survival (can be changed)
#' @param rna_data,clinical_data,what_to_check,gene_identifier
#' @keywords data_processing
#' @export
#' @examples
#' prepareRNAdata()

prepareRNAdata <- function(rna_data,clinical_data,what_to_check,gene_identifier){
  rna_data <- na.omit(rna_data)
  if (gene_identifier == 'Entrez'){
    rna_data$Hugo_Symbol <- NULL
  }
  else if (gene_identifier == 'Hugo'){
    rna_data$Entrez_Gene_Id<- NULL
  }
  rownames(rna_data) <- make.names(rna_data[,1], unique=TRUE)
  rna_data <- rna_data[,-1]
  rna_data2 = t(rna_data)
  dummy_col <- matrix(0,length(rna_data2[,1]),1)
  Gene_names <- as.character(rownames(rna_data2))
  index_vector <- c()
  for (i in 1:length(clinical_data$PATIENT_ID)){
    idx <- which(as.factor(rownames(rna_data2)) == clinical_data$PATIENT_ID[i])
    if (length(idx) != 0){
      index_vector <- c(index_vector,idx)
    }
  }
  #id_clinical = which(as.factor(rownames(rna_data2)) %in% clinical_data$PATIENT_ID) # Could have to re-evaluate this...
  id_clinical = index_vector
  dummy_col2 <- matrix(NA,length(rna_data2[,1]),1)
  if (what_to_check == "survival"){
    print("Going for the deaths....")
    dummy_col2[id_clinical] = as.character(clinical_data$OS_MONTHS[id_clinical])
    Survival_months = suppressWarnings(as.numeric(as.character(dummy_col2)))
    three_groups <- as.factor(as.numeric(cut_number(Survival_months,3)))
    rna_data2 = cbind.data.frame(Survival_months,rna_data2)
    rna_data2 = cbind.data.frame(three_groups,rna_data2)
    rna_data2 = cbind.data.frame(Gene_names,rna_data2)
    rna_data2 <- na.omit(rna_data2)
    #rna_data2 <- rna_data2[ - as.numeric(which(apply(rna_data2, 2, var) == 0))] # remove zero variance columns
    levels(rna_data2$three_groups) <- list('Low Survival'=1, 'Moderate Survival'=2, 'High Survival'=3)
  }
  else if (what_to_check == "race"){
    print("Going for the races....")
    dummy_col2[1:length(clinical_data$RACE[id_clinical])] = as.character(clinical_data$RACE[id_clinical])
    races = suppressWarnings(as.factor(as.character(dummy_col2)))
    rna_data2 = cbind.data.frame(races,rna_data2)
    rna_data2 = cbind.data.frame(Gene_names,rna_data2)
    rna_data2 <- na.omit(rna_data2)
    #rna_data2 <- rna_data2[ - as.numeric(which(apply(rna_data2, 2, var) == 0))] # remove zero variance columns
    rna_data2 = subset(rna_data2, races!="[Not Available]")
    rna_data2 = subset(rna_data2, races!="C50.9")
    rna_data2 = subset(rna_data2, races!="DiseaseFree")
  }
  return(rna_data2)
}
