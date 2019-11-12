#' # A function for extracting the top/bottom n% of a gene expression matrix based on the expression of 'list_of_genes'
#'
#' @param dataset,genes,nth_percentage,cor_data2,type_of_effect,hazard_or_pval,print_plot
#' @keywords data_processing
#' @export
#' @examples
#' segregate_dataset()

segregate_dataset <- function(dataset,genes,nth_percentage,cor_data2,type_of_effect,hazard_or_pval,print_plot=TRUE){
  Packages <- c("survminer", "survival")
  lapply(Packages, require, character.only = TRUE)
  if (missing(print_plot))
  {
    print_plot <- TRUE
  }
  if (missing(hazard_or_pval))
  {
    hazard_or_pval <- 'pval'
  }
  if (print_plot == TRUE){
    id_no_effect <- which(type_of_effect == 'No effect!')
  } else {
    id_no_effect <- c()
  }
  if (class(genes) != 'data.frame'){
    genes <- data.frame(genes)
    colnames(genes) <- 'Group1'
    if(length(id_no_effect) != 0){
      type_of_effect <- type_of_effect[-id_no_effect]
      genes <- data.frame(genes[-id_no_effect,])
      colnames(genes) <- 'Group1'
    } else {
      type_of_effect <- type_of_effect
      genes <- genes
    }
  } else {
    if(length(id_no_effect) != 0){
      type_of_effect <- type_of_effect[-id_no_effect]
      genes <- genes[-id_no_effect,]
    } else {
      type_of_effect <- type_of_effect
      genes <- genes
    }
  }
  #print(genes[,1])
  id_used_3 <- which(rownames(dataset) %in% genes[,1])
  if (length(id_used_3) == 0){
    print('No gene with any significant effect on survival, this will result in an error..')
  }
  #print(id_used_3)
  dataset2 <- dataset[id_used_3,]
  #print(sd(as.numeric(dataset2[1096,])))
  #print(dataset2)
  zdata <- standardize(dataset2)
  # Add a step where the NA values (I dont know why they actually exist) are turned into 0.
  # ADDED 090519 as a test-------------------
  for (i in 1:dim(dataset2)[1]){
      id_gene <- which(genes$Group1 %in% rownames(dataset2)[i])
      effect_type <- type_of_effect[id_gene]
      if (length(effect_type == 'Beneficial!') == 0){
        zdata[i,] = zdata[i,]
        print('No change in scoring...')
      } else if (length(effect_type == 'Beneficial!') == 1){
        if (effect_type == 'Beneficial!'){
          zdata[i,] = -as.numeric(as.character(zdata[i,]))
        } else {
          zdata[i,] = zdata[i,]
        }
      }
  }
  if (all(is.na(c(NA, zdata)))){ # test if all values are the same (no sd = NA values in standardized data)
    if (print_plot == T){
      p3 <- 'Standard deviation is 0, no standardization possible!'
    } else if (print_plot == F){
      p3 <- 1
    }
  } else {
    zdata[is.na(zdata)] <- 0
    sum_genes = colSums(zdata) # This could use some adjustment, 1-2 genes could be carrying the rest.
    #print(sum_genes)
    thresh_val_top = quantile(sum_genes,1-nth_percentage)
    patients_top <- data.frame(names(sum_genes)[which(sum_genes >= thresh_val_top)])
    colnames(patients_top) = 'patients'
    thresh_val_bottom = quantile(sum_genes,nth_percentage)
    patients_bottom <- data.frame(names(sum_genes)[which(sum_genes <= thresh_val_bottom)])
    colnames(patients_bottom) = 'patients'
    dataset_top <- dataset[,which(colnames(dataset) %in% patients_top$patients)]
    dataset_bottom <- dataset[,which(colnames(dataset) %in% patients_bottom$patients)]
    # Here I need to assert whether cor_data2 contains a 'patient','survival' and 'status' vectors. Otherwise, it will
    # result in a boring error that is hard to track unless going directly through all lines in the function...
    df_patients_top <- cor_data2[which(cor_data2$patient %in% colnames(dataset_top)),]
    df_patients_top <- cbind(df_patients_top,rep('top',nrow(df_patients_top)))
    colnames(df_patients_top)[ncol(df_patients_top)] = 'Expression'
    df_patients_bottom <- cor_data2[which(cor_data2$patient %in% colnames(dataset_bottom)),]
    df_patients_bottom <- cbind(df_patients_bottom,rep('bottom',nrow(df_patients_bottom)))
    colnames(df_patients_bottom)[ncol(df_patients_bottom)] = 'Expression'
    df_fin = rbind(df_patients_top,df_patients_bottom)
    # Now do the survival curve:
    fit_chosen <- survfit(Surv(survival, status) ~ Expression, data = df_fin)
    
    if (print_plot == TRUE){
      if (dim(dataset2)[1] == 0) { # Checking if the dataframe is empty (if checking one gene's effect on survival)
        print('No genes left in dataframe!')
      }
      p3 <- ggsurvplot(fit_chosen, data = df_fin, risk.table=T,pval=T, conf.int=T,
                       xlim = c(0,60),xlab = 'Time in months', # Choose length of follow-up time (120 = 10 years, 60 = 5 years)
                       break.time.by = 12,
                       legend.labs = c('Favourable signature','Unfavourable signature'),
                       #palette = c("gray72", "gray24"),
                       palette = c('gray10','darkred'),
                       ggtheme = theme_classic2())
      p3 <- ggpar(p3,font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"))
      if (length(which(colnames(df_fin) == 'age')== 1) & length(which(colnames(df_fin) == 'stage')== 1)){
        df_fin$stage <- as.numeric(df_fin$stage)
        covariates <- c("Expression","age","stage")
      } else {
        covariates <- c("Expression")
      }
      univ_formulas <- sapply(covariates,
                              function(x) as.formula(paste('Surv(survival, status)~', x)))
      univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df_fin)})
      univ_results <- lapply(univ_models,
                             function(x){ 
                               x <- summary(x)
                               p.value<-signif(x$wald["pvalue"], digits=2)
                               #p.value <- x$sctest[3]
                               wald.test<-signif(x$wald["test"], digits=2)
                               beta<-signif(x$coef[1], digits=2);#coeficient beta
                               HR <-signif(x$coef[2], digits=2);#exp(beta)
                               HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                               HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                               HR <- paste0(HR, " (", 
                                            HR.confint.lower, "-", HR.confint.upper, ")")
                               res<-c(beta, HR, wald.test, p.value)
                               names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                             "p.value")
                               return(res)
                             })
      res <- t(as.data.frame(univ_results, check.names = FALSE))
      print(as.data.frame(res))
    } else {
      if (dim(dataset2)[1] == 0) { # Checking if the dataframe is empty (if checking one gene's effect on survival)
        print('No genes left in dataframe!')
        p3 <- NA
      } else {
        res.cox_cand <- coxph(Surv(survival, status) ~ Expression, data = df_fin)
        hehe <- summary(res.cox_cand)
        if (length(which(colnames(df_fin) == 'age')== 1) & length(which(colnames(df_fin) == 'stage')== 1)){
          df_fin$stage <- as.numeric(df_fin$stage)
          covariates <- c("Expression","age","stage")
        } else {
          covariates <- c("Expression")
        }
        univ_formulas <- sapply(covariates,
                                function(x) as.formula(paste('Surv(survival, status)~', x)))
        univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df_fin)})
        univ_results <- lapply(univ_models,
                               function(x){ 
                                 x <- summary(x)
                                 p.value<-signif(x$wald["pvalue"], digits=2)
                                 #p.value <- x$sctest[3]
                                 wald.test<-signif(x$wald["test"], digits=2)
                                 beta<-signif(x$coef[1], digits=2);#coeficient beta
                                 HR <-signif(x$coef[2], digits=2);#exp(beta)
                                 HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                 HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                                 HR2 <- paste0(HR, " (", 
                                               HR.confint.lower, "-", HR.confint.upper, ")")
                                 res<-c(beta, HR, HR2, wald.test, p.value)
                                 names(res)<-c("beta","HR","95% CI for HR", "wald.test", 
                                               "p.value")
                                 return(res)
                               })
        res <- t(as.data.frame(univ_results, check.names = FALSE))
        #print(as.data.frame(res))
        #print(res)
        #p3 <- hehe$sctest[3] # This is the p-value - return that!
        # For the analysis of 20.09.19 - Return hazard ratios instead:
        if (hazard_or_pval == 'pval'){
          p3 <- hehe$sctest[3] # p-value from the log-rank test...
        } else {
          p3 <- as.numeric(res[1,2]) # the HR...
        }
      }
    }
  }
  return(p3)
}