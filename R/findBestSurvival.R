#' Extracting the best split of survival of patients from expression of genes.
#'
#' A function for optimizing the value chosen for survival signature generation
#' @param dataset,genelist,quantile_list,cor_data,print_plot,effect_type
#' @keywords survival
#' @export
#' @examples
#' findBestSurvival()

findBestSurvival <- function(dataset,genelist,quantile_list,cor_data,print_plot=TRUE,effect_type=TRUE){
  if (colnames(genelist) != 'Group1'){
    colnames(genelist) = 'Group1'
  }
  if (missing(effect_type))
  {
    effect_type <- FALSE
  }
  if (missing(print_plot))
  {
    print_plot <- FALSE
  } else if (effect_type == TRUE){
    print_plot <- FALSE
  }
  p_comp = Inf # start with the highest...
  lowest_quantile = NULL # Start with the lowest quantile
  q_vector <- quantile_list
  cor_data2 <- cor_data
  for (i in 1:length(q_vector)){
    patients_assignment_vector <- prepareSurvivalDataTCGA(dataset,genelist,q_vector[i],0.75)
    cor_data4 <- cbind(cor_data2,as.factor(patients_assignment_vector))
    colnames(cor_data4)[ncol(cor_data4)] <- c('Cl3_high_low')
    res.cox2 <- coxph(Surv(survival, status) ~ Cl3_high_low, data = cor_data4)
    sum_coxph <- summary(res.cox2)
    pval_to_test <- sum_coxph$sctest[3]
    print(pval_to_test)
    if (pval_to_test < p_comp){
      p_comp <- pval_to_test
      lowest_quantile <- q_vector[i]
      cor_data_chosen <- cor_data4
    } else if (pval_to_test >= p_comp) {
      p_comp <- p_comp
      lowest_quantile <- lowest_quantile
      cor_data_chosen <- cor_data_chosen
      lowest_quantile <- lowest_quantile
    }
  }
  # Now output the correct graph:
  print(lowest_quantile)
  patients_assignment_vector <- prepareSurvivalDataTCGA(dataset,genelist,lowest_quantile,0.75)
  cor_data_chosen <- cbind(cor_data2,as.factor(patients_assignment_vector))
  colnames(cor_data_chosen)[ncol(cor_data_chosen)] <- c('Cl3_high_low')
  res.cox_cand <- coxph(Surv(survival, status) ~ Cl3_high_low, data = cor_data_chosen)
  type_of_effect <- sign(res.cox_cand$coefficients)
  fit_chosen <- survfit(Surv(survival, status) ~ Cl3_high_low, data = cor_data_chosen)
  if(print_plot == TRUE){
    p3 <- ggsurvplot(fit_chosen, data = cor_data_chosen, risk.table=T,pval=T, conf.int=T,
                     xlim = c(0,120),xlab = 'Time in months',
                     break.time.by = 12,
                     legend.labs = c('Higher expression of genes','Lower expression of genes'),
                     ggtheme = theme_light())
    p3 <- ggpar(p3,font.x = c(14, "bold", "black"),
                font.y = c(14, "bold", "black"))
    print(paste0("Q-value chosen is ", lowest_quantile))
  } else if (print_plot == FALSE){
    if (effect_type == FALSE){
      hehe <- summary(res.cox_cand)
      p3 <- hehe$sctest[3] # This is the p-value - return that!
    } else {
      if (type_of_effect == 1){
        p3 <- 'Beneficial!'
      } else if (type_of_effect == -1){
        p3 <- 'Harmful!'
      }
    }
  }
  return(p3)
}
