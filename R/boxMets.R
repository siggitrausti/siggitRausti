#' Create beautiful boxplots from data. 
#'
#' @param dataset,id_mets,x_variable,fill_variable,common_legend
#' @keywords data_processing
#' @export
#' @examples
#' boxMets()

boxMets <- function(dataset,id_mets,x_variable,fill_variable,common_legend=T){
  # SINCE I HAVE NOT ADDED THE 'PLOT_METABOLITES' FUNCTION TO SIGGITRAUSTI PACKAGE, THIS SHORTCUT WILL HAVE TO DO..
  Packages <- c("ggpubr", "gplots","cowplot","gridExtra","siggitRausti")
  lapply(Packages, require, character.only = TRUE)
  if (missing(id_mets))
  {
    for (k in 1:dim(dataset)[2]){
      if (!is.numeric(dataset[,k])){
        k = k+1
      } else {
        id_mets <- k:dim(dataset)[2]
      }
    }
  }
  if (missing(common_legend))
  {
    common_legend=T
  }
  plot_list <- list()
  for (i in 1:length(id_mets)){
    p1 <- ggboxplot(dataset, x = x_variable, y = paste0("`", colnames(dataset)[id_mets[i]], "`"),
                    fill = fill_variable, palette = ColBrew('JCO'),
                    ylab = "Metabolite", xlab = "Treatment",
                    width = 0.7,size=0.9)
    p1 <- plotLookForPaper(p1,'Ratios','',rotate_check=T,legend_check = F)
    p1 <- p1 + ggtitle(paste(colnames(dataset)[id_mets[i]]))
    plot_list[[i]] <- p1
  }
  p_fin <- plot_metabolites(plot_list,nplots = length(id_mets),common_legend = common_legend)
  return(p_fin)
}
  