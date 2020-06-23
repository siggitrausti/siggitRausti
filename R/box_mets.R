#' Create beautiful boxplots from data. 
#'
#' @param dataset,id_mets,x_variable,fill_variable,common_legend,stat_test,ref_group
#' @keywords data_processing
#' @export
#' @examples
#' box_mets()

box_mets <- function(dataset,id_mets,x_variable,fill_variable,common_legend=T,stat_test=F,ref_group,y_title){
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
  if (missing(stat_test))
  {
    stat_test=F
  }
  if (missing(ref_group))
  {
    ref_group='.all.'
  }
  if (missing(y_title)){
    y_title = 'Ratios'
  }
  plot_list <- list()
  for (i in 1:length(id_mets)){
    if (stat_test == T){
      my_ref_group <- ref_group
      p1 <- ggboxplot(dataset, x = x_variable, y = paste0("`", colnames(dataset)[id_mets[i]], "`"),
                      fill = fill_variable, palette = ColBrew('MetaboA'),
                      ylab = "Metabolite", xlab = "Treatment",
                      width = 0.7,size=0.9) + 
        stat_compare_means(size=10,label='p.signif',hide.ns = T,ref.group = my_ref_group)
      p1 <- plotLookForPaper(p1,y_title,'',rotate_check=T,legend_check = F)
      p1 <- p1 + ggtitle(paste(colnames(dataset)[id_mets[i]]))
      p1 <- p1 + theme(axis.text.x = element_text(face="bold", color="grey10", 
                                                  size=18),
                       axis.text.y = element_text(face="bold", color="grey10", 
                                                  size=14),
                       axis.title.x = element_text(face="bold", color="grey10", 
                                                   size=18),
                       axis.title.y = element_text(face="bold", color="grey10", 
                                                   size=18))
      p1 <- p1 + theme(panel.grid.major.y = element_line(color='grey'))
      plot_list[[i]] <- p1
    } else {
      p1 <- ggboxplot(dataset, x = x_variable, y = paste0("`", colnames(dataset)[id_mets[i]], "`"),
                      fill = fill_variable, palette = ColBrew('Klimt'),
                      ylab = "Metabolite", xlab = "Treatment",
                      width = 0.7,size=0.9)
      p1 <- plotLookForPaper(p1,y_title,'',rotate_check=T,legend_check = F)
      p1 <- p1 + ggtitle(paste(colnames(dataset)[id_mets[i]]))
      p1 <- p1 + theme(axis.text.x = element_text(face="bold", color="grey10", 
                                                  size=10),
                       axis.text.y = element_text(face="bold", color="grey10", 
                                                  size=10),
                       axis.title.x = element_text(face="bold", color="black", 
                                                   size=14),
                       axis.title.y = element_text(face="bold", color="black", 
                                                   size=14))
      p1 <- p1 + theme(panel.grid.major.y = element_line(color='grey'))
      plot_list[[i]] <- p1
    }
    
  }
  p_fin <- plot_metabolites(plot_list,nplots = length(id_mets),common_legend = common_legend)
  return(p_fin)
}

  