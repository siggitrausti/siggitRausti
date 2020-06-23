#' Create beautiful boxplots from data. 
#'
#' @param houston_test,return_val
#' @keywords plots
#' @export
#' @examples
#' ROC_maker()

ROC_maker <- function(houston_test,return_val){
  # A function to make a ROC classifier from a dataset. All variables are taken into account for prediction of binary response. 
  # The binary response should be the first column, and then the rest should be the variables....
  if (missing(return_val)){
    return_val = 'plot'
  }
  
  packages_needed = c('mlr','pROC','stringr','ggplot2')
  lapply(packages_needed, require, character.only = TRUE)
  target_variable = colnames(houston_test)[1]
  task <- makeClassifTask(data = houston_test, target = target_variable)
  glm.learner <- mlr::makeLearner("classif.binomial",predict.type = 'prob')
  #logistic.learner <- makeLearner("classif.binomial",predict.type = "prob")
  fit <- mlr::train(glm.learner,task)
  pred <- predict(fit, task)
  #names_for_legend = colnames(houston_test)[2:ncol(houston_test)]
  #names_for_legend <- str_sub(names_for_legend, 2, -2)
  roc <- generateThreshVsPerfData(pred, list(fpr, tpr))
  #plotROCCurves(roc) + theme_minimal() + coord_equal()
  obj <- pROC::roc(pred$data$truth, pred$data[[3]], ci=TRUE, plot=FALSE)
  print(paste0('AUC of model: ',signif(obj$auc,2), ' (95% CI ',signif(obj$ci[1],2),' to ',signif(obj$ci[3],2),')'))
  #print(paste0('Sensitivity of model: ',signif(pred$data$tpr)))
  roc_with_ci <- function(obj) {
    ciobj <- ci.se(obj, specificities = seq(0, 1, l = 25))
    dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                         lower = ciobj[, 1],
                         upper = ciobj[, 3])
    
    ggroc(obj,size=1) +
      theme_minimal() +
      geom_abline(
        slope = 1,
        intercept = 1,
        linetype = "dashed",
        alpha = 0.7,
        color = "grey10",
        size = 1) + 
      coord_equal() +
      geom_ribbon(
        data = dat.ci,
        aes(x = x, ymin = lower, ymax = upper),
        fill = "steelblue",
        alpha = 0.2) + 
      #ggtitle(capture.output(obj$ci)) + 
      ylab('Sensitivity') + xlab('Specificity') + 
      theme(plot.title   = element_text(size=20, hjust= .5, vjust = 2, face = "bold"),
            axis.text.x = element_text(face="bold", color="grey10", 
                                       size=12),
            axis.text.y = element_text(face="bold", color="grey10", 
                                       size=12),
            axis.title.x = element_text(face="bold", color="grey10", 
                                        size=18),
            axis.title.y = element_text(face="bold", color="grey10", 
                                        size=18), 
            #legend.title = element_blank(),
            #legend.spacing.y = unit(0, "mm"),
            panel.border = element_rect(colour = "white", fill=NA),
            aspect.ratio = 0.9, axis.text = element_text(colour = 1, size = 12),
            #legend.background = element_blank(),
            #legend.box.background = element_rect(colour = "black"),
            axis.line = element_line(colour ="black", size = 0.8),
            axis.ticks = element_line(colour ="black", size = 0.8),
            legend.text  = element_text(size =14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"))
  } 
  if (return_val == 'AUC'){
    return(obj$auc)
  } else if (return_val == 'plot'){
    my_plot <- roc_with_ci(obj)
    my_plot <- my_plot + annotate('text',x=0.3,y = 0.2,label=paste0('AUC of model: ',signif(obj$auc,2), '\n (95% CI ',signif(obj$ci[1],2),' to ',signif(obj$ci[3],2),')'),
                                  size = 6)
    
    return(my_plot)
  }
}