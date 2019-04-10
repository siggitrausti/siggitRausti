#' Color brewer
#'
#' This function allows you to generate a color palette for images
#' @param type,ncolors
#' @keywords colors
#' @export
#' @examples
#' ColBrew()
#'
ColBrew <- function(type, ncolors){

  outp = NULL;
  nargin <- length(as.list(match.call())) -1;
  if (missing(ncolors)){
    ncolors <- 9
  }

 if (ncolors > 9) {
    print("No more than 9 colors can be selected.")
    return(outp)
  }

  PI      <- c("#1e434c", "#CB3618", "#2294D5","white","white",
               "white","white","white","white")

  Google  <- c("#0057e7","#d62d20","#008744","#ffa700","white",
                "white", "white", "white", "white")

  Klimt   <- c("#0B3C5D","#D9B310","#328CC1","#7A3200", "#436229",
               "#DA888C", "#CB1E00", "#E29D50", "black")

  Gogh    <- c("#527192","#F0BF00","#007849","#262228", "#6D7255",
               "#974F36","#F0BF00","#6B8546","#CBAB5E")

  Warhol  <- c("#ef0000","#3FDBE6","#ffdb00","#62ce72",'#7AA6DC',
               "white", "white", "white", "white")

  JCO     <- c('#0073C2', '#CD534C','#EFC000', '#868686','#003C67',
               '#8F7700', '#3B3B3B', '#A73030', '#4A6990')

  NEJM     <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF", "#7876B1FF",
                "#6F99ADFF", "#FFDC91FF","#EE4C97FF","white")

  startrek <- c("#df0000", "#0099f6", "#f2c300", "#00b844","#FFCD00FF",
                "#7C878EFF","#00B5E2FF","#00AF66FF","white")

  grey     <- c("#000000","#a7a7a7",'#575757',"#d1d1d12","#666666",
                "#333333","#999999","#111111","white")

  NPG      <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF",'#F39B7F',
                '#8491B4','#91D1C2', '#DC0000', '#7E6148')

  AAAS     <- c("#3B4992FF","#EE0000FF","#008B45FF","#631879FF",'#5F559B',
                '#BB0021', '#008280', '#A20056', '#808180')

  Lancet   <- c('#00468B', '#ED0000', '#42B540', '#0099B4', '#925E9F',
                '#FDAF91', '#AD002A', '#ADB6B6', '#1B1919')

  JAMA     <- c("#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF", "#79AF97FF",
                "#6A6599FF", "#80796BFF","white","white")

  D3       <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF")

  Simpsons <- c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF",
                "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#71D0F5FF")

  Futurama <- c("#FF6F00FF", "#C71000FF", "#008EA0FF", "#8A4198FF", "#5A9599FF",
                "#FF6348FF", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF")

  Schwifty <- c("#FAFD7CFF", "#82491EFF", "#B7E4F9FF", "#24326FFF", "#FB6467FF",
                "#526E2DFF", "#E762D7FF", "#E89242FF", "#FAE48BFF")

  Mononoke <- c("#4870B6", "#DD6143", "#A69A69", "#3E3637", "#9BAFD1",
                "#6A7A5A", "#ECCE94", "#F3FAFC", "#FAE48BFF")

  Toystory <- c("#0092E7", "#EC4D4A", "#F9D401", "#A3BB4B", "#9A5E9A",
                "#FAC893", "#F9CCE3", "#A84327", "#F98B5A")

  Snowwhite <- c("#1F3984", "#E61616", "#F2C148", "#5A997C", "#7BA8C9",
                "#C06B66", "#F38F3A", "#F9E58F", "#EDEFF4")
  MetaboA   <- c("red", "green3", "blue", "cyan", "magenta", "yellow","gray","black","firebrick3")


  ColorSets <- data.frame(PI,Google, Klimt, Gogh, Warhol, JCO, NEJM, startrek, grey, NPG, AAAS,Lancet, JAMA,
                          D3, Simpsons, Futurama, Schwifty, Mononoke, Toystory, Snowwhite, MetaboA);

  colID <- which(colnames(ColorSets) == type);

  outp <- ColorSets[1:ncolors,colID]
  outp <- as.character(outp)
  return(outp)
}
