#' A function to plot multiple plots from a list to the same graph, with a common legend or not
#'
#' @param plots_metabolites,nplots,common_legend
#' @keywords data_processing
#' @export
#' @examples
#' plot_metabolites()

plot_metabolites <- function(plots_metabolites,nplots,common_legend = T) {
 # A function to plot multiple plots from a list to the same graph, with a common legend or not 
  # Basicly works for length of metabolites from 4-15, can be changed to work for other numbers...
  require('gridExtra','ggpubr','ggplot2')
  if (missing(common_legend))
  {
    common_legend <- FALSE
  }
  if (common_legend == T){
    if (nplots == 1){
      outp <- ggarrange(plots_metabolites[[1]],
                        ncol = 1, nrow = 1,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 2){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]],
                        ncol = 2, nrow = 1,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 3){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]],plots_metabolites[[3]],
                        ncol = 2, nrow = 2,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 4){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]], plots_metabolites[[4]],
                        ncol = 2, nrow = 2,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 5){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]],
                        ncol = 3, nrow = 2,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 6){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]],
                        ncol = 3, nrow = 2,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 7){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]],
                        ncol = 3, nrow = 3,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 8){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        ncol = 3, nrow = 3,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 9){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]],
                        ncol = 3, nrow = 3,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 10){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]], 
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]],
                        ncol = 3, nrow = 4,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 11){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]],
                        ncol = 3, nrow = 4,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 12){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        ncol = 3, nrow = 4,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 13){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        plots_metabolites[[13]], 
                        ncol = 3, nrow = 5,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 14){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        plots_metabolites[[13]], plots_metabolites[[14]], 
                        ncol = 3, nrow = 5,
                        common.legend = TRUE, legend = "bottom")
    } else if (nplots == 15){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        plots_metabolites[[13]], plots_metabolites[[14]], plots_metabolites[[15]], 
                        ncol = 3, nrow = 5,
                        common.legend = TRUE, legend = "bottom")
    }
  } else if (common_legend == F){
    if (nplots == 1){
      outp <- ggarrange(plots_metabolites[[1]],
                        ncol = 1, nrow = 1,
                        common.legend = FALSE)
    } else if (nplots == 2){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]],
                        ncol = 2, nrow = 1,
                        common.legend = FALSE)
    } else if (nplots == 3){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]],plots_metabolites[[3]],
                        ncol = 2, nrow = 2,
                        common.legend = FALSE)
    } else if (nplots == 4){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]], plots_metabolites[[4]],
                        ncol = 2, nrow = 2,common.legend = FALSE)
    } else if (nplots == 5){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]],
                        ncol = 3, nrow = 2,common.legend = FALSE)
    } else if (nplots == 6){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]],
                        ncol = 3, nrow = 2,common.legend = FALSE)
    } else if (nplots == 7){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]],
                        ncol = 3, nrow = 3,common.legend = FALSE)
    } else if (nplots == 8){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        ncol = 3, nrow = 3,common.legend = FALSE)
    } else if (nplots == 9){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]],
                        ncol = 3, nrow = 3,common.legend = FALSE)
    } else if (nplots == 10){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]], 
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]],
                        ncol = 3, nrow = 4,common.legend = FALSE)
    } else if (nplots == 11){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]],
                        ncol = 3, nrow = 4,common.legend = FALSE)
    } else if (nplots == 12){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        ncol = 3, nrow = 4,common.legend = FALSE)
    } else if (nplots == 13){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        plots_metabolites[[13]], 
                        ncol = 3, nrow = 5,common.legend = FALSE)
    } else if (nplots == 14){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        plots_metabolites[[13]], plots_metabolites[[14]], 
                        ncol = 3, nrow = 5,common.legend = FALSE)
    } else if (nplots == 15){
      outp <- ggarrange(plots_metabolites[[1]], plots_metabolites[[2]], plots_metabolites[[3]],plots_metabolites[[4]],
                        plots_metabolites[[5]], plots_metabolites[[6]], plots_metabolites[[7]], plots_metabolites[[8]], 
                        plots_metabolites[[9]], plots_metabolites[[10]], plots_metabolites[[11]], plots_metabolites[[12]], 
                        plots_metabolites[[13]], plots_metabolites[[14]], plots_metabolites[[15]], 
                        ncol = 3, nrow = 5,common.legend = FALSE)
    }
  }
  return(outp)
}