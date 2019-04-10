#' Making plots look good
#'
#' @param p,y_text,x_text,rotate_check,legend_check
#' @keywords plots
#' @export
#' @examples
#' plotLookForPaper()


plotLookForPaper <- function(p,y_text,x_text,rotate_check=FALSE,legend_check=TRUE){
  if (missing(rotate_check))
  {
    rotate_check <- FALSE
  }
  if (missing(legend_check))
  {
    legend_check <- TRUE
  }
  if (missing(x_text))
  {
    x_text = ''
  }
  p <- p + theme_bw()
  p <- p + theme(plot.title   = element_text(size=20, hjust= .5, vjust = 2, face = "bold"),
                 legend.title = element_blank(),
                 legend.spacing.y = unit(0, "mm"),
                 panel.border = element_rect(colour = "white", fill=NA),
                 aspect.ratio = 0.9, axis.text = element_text(colour = 1, size = 12),
                 legend.background = element_blank(),
                 legend.box.background = element_rect(colour = "black"),
                 axis.title.y = element_text(size=16, face = "bold", colour = "black",
                                             margin = margin(t = 0, r = 15, b = 0, l = 0)),
                 axis.title.x = element_text(size=16, face = "bold", colour = "black",
                                             margin = margin(t = 15, r = 0, b = 15, l = 0)),
                 axis.text.y  = element_text(size=14, face = "bold", colour = "black"),
                 axis.line = element_line(colour ="black", size = 1),
                 axis.ticks = element_line(colour ="black", size = 1),
                 legend.text  = element_text(size =14),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  p <- p + ylab(y_text) + xlab(x_text)
  if(rotate_check==TRUE)
  {
    p <- p + theme(axis.text.x  = element_text(size=14, face = "bold", colour = "black",angle=45, hjust=1))
    }
  else {
    p <- p + theme(axis.text.x  = element_text(size=14, face = "bold", colour = "black"))
    }
  if(legend_check==FALSE)
  {
    p <- p + theme(legend.position  = "none")
  }
  return(p)
}
