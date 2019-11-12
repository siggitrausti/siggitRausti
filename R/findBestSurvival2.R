findBestSurvival2 <- function(dataset,genelist,quantile_list,cor_data,print_plot=TRUE,effect_type=TRUE){
  if (class(genelist) != 'data.frame'){
    genelist <- data.frame(genelist)
  }
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
  if (length(q_vector) > 1){
    for (i in 1:length(q_vector)){
      patients_assignment_vector <- prepareSurvivalDataTCGA(dataset,genelist,q_vector[i],0.6) # This 0.75 value needs optimization!!
      if (all(patients_assignment_vector == 0)){
        lowest_quantile = q_vector[i]
        break
      } else {
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
    }
  } else {
    lowest_quantile = q_vector
  }
  # Now output the correct graph:
  patients_assignment_vector <- prepareSurvivalDataTCGA(dataset,genelist,lowest_quantile,0.6)
  cor_data_chosen <- cbind(cor_data2,as.factor(patients_assignment_vector))
  colnames(cor_data_chosen)[ncol(cor_data_chosen)] <- c('Cl3_high_low')
  if (all(patients_assignment_vector == 0)){
    type_of_effect <- 0
    if (print_plot == TRUE){
      p3 <- 'No plot will be plotted, since there is no significant segregation'
      print(p3)
    } else if (print_plot == FALSE){
      if (effect_type == FALSE){
        p3 <- 1
      } else {
        if (type_of_effect == 0){
          p3 <- 'No effect!'
        }
      }
    }
  } else {
    res.cox_cand <- coxph(Surv(survival, status) ~ Cl3_high_low, data = cor_data_chosen)
    type_of_effect <- sign(res.cox_cand$coefficients)
    fit_chosen <- survfit(Surv(survival, status) ~ Cl3_high_low, data = cor_data_chosen)
    if(print_plot == TRUE){
      p3 <- ggsurvplot(fit_chosen, data = cor_data_chosen, risk.table=T,pval=T, conf.int=F,
                       xlim = c(0,120),xlab = 'Time in months',
                       break.time.by = 12,
                       legend.labs = c('High expression','Low expression'),
                       palette = c('darkred','gray10'),
                       ggtheme = theme_classic2())
      p3 <- ggpar(p3,font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"))
      print(paste0("Q-value chosen is ", lowest_quantile))
    } else if (print_plot == FALSE){
      if (effect_type == FALSE){
        hehe <- summary(res.cox_cand)
        p3 <- hehe$sctest[3] # This is the p-value - return that!
      } else {
        hehe <- summary(res.cox_cand)
        p_value2 <- hehe$sctest[3]
        if (type_of_effect == 1 & p_value2 < 0.05){
          p3 <- 'Beneficial!'
        } else if (type_of_effect == -1 & p_value2 < 0.05){
          p3 <- 'Harmful!'
        } else {
          p3 <- 'No effect!'
        }
      }
    }
  }
  return(p3)
}
