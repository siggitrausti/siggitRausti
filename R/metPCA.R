#' PCA plotting of data:
#'
#' metPCA from Freyr Johannson.
#' @param dataset,groups,scale,pcs,center,log,power,colors,CI,label,corplot,minvalue,narm,scree,nh,pointsize,corrList,loadinglist
#' @keywords PCA
#' @export
#' @examples
#' metPCA()

metPCA <- function(dataset, groups, scale, pcs ,center, log, power,colors, CI,label, corplot, minvalue, narm ,scree, nh,
                   pointsize, corrList,loadinglist)
{


# Just to make sure everything is in order regarding libraries:
packages = c("ggplot2","scales", "ggrepel") # Removed  "MetabolAnalyze"

  #use this function to check if each package is on the local machine
  #if a package is installed, it will be loaded
  #if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })

output <- NULL

# Should possibly remove....
if (!exists("ColBrew"))
{
   if (!file.exists("ColBrew.R")) {
      cat("You need the ColBrew function, talk to Freyr or modify this function.\n")
      cat("Maybe it's in the wrong folder.")
      return()
   }
   else {
     source("ColBrew.R")
   }
}

# Check the optional parametes:
if (missing (colors))
  {
    colors = ColBrew('Futurama',9)
  }

if (missing(log))
  {
    log = FALSE
  }
if (missing(power))
  {
  power = FALSE
  }
if (missing(scale))
  {
   scale = c("none")
  }
if (missing(center))
  {
   center = TRUE
  }

if(missing(nh))
  {
   nh = 10
  }

bLabel <- TRUE
if (missing(label))
  {
   bLabel <- FALSE
  }

if (missing(corplot))
  {
  corplot <- FALSE
  }

if (missing(pointsize))
  {
  pointsize = 3.5
}
if (missing(pcs))
  {
    pcs = c(1,2)
  }
if (missing(groups))
{
  groups = c("none")
}
if (missing(scree))
{
  scree = F
}
if (missing(CI)){
  CI = FALSE
}
if (missing(loadinglist))
{
  loadinglist = FALSE
}
if (missing(corrList))
{
  corrList = FALSE
}
if (missing(minvalue))
{
  minvalue = FALSE
}
if (missing(narm))
{
  narm = T
}


# Scaling functions, no need to do auto since scale is a base function in R

pareto <- function(dataf)
{
  dataf <- apply(dataf,2, function(x) (x-mean(x))/sqrt(sd(x)))
  return(dataf)
}

Range <- function(dataf)
{
  dataf <- apply(dataf,2, function(x) (x-mean(x)) / (max(x,na.rm = T)-min(x,na.rm = T)))

  return(dataf)
}

Vast <- function(dataf)
{
  dataf <- apply(dataf,2, function(x) (x-mean(x))/sd(x) * mean(x)/sd(x))
  return(dataf)
}

LevelScaling <- function(dataf)
{
  dataf <- apply(dataf,2,function(x) (x-mean(x))/mean(x) )
  return(dataf)
}

# The standard form of metabolic data has the first x columns as categorical
# with time, or another numerical column not used in the analysis as the
# first numerical column, every column after that should be numerical and
# will be included in the analysis (check out loadData.r)

# Find the first numeric column (normally the "Time" column)
i = 1;
while (class(dataset[,i]) != "numeric")
{
  i = i+1;
}

pcadata <- dataset[,-(1:i)] # record the numerical values used in the pca
catdata <- dataset[,(1:i)]  # record the categorical values.

# Replace missing values with the half minimum value in the dataset,
# Which is the default method in metaboanalyst.

# Lets get rid of zeros and negative numbers at the same time:
# pcadata[pcadata <= 0] <- NA
# There are possibly still negative numbers, do not want to use them as
# the mininum value.
if (minvalue)
{
  pcadata[pcadata <= 0] <- NA
  minValue <- min(pcadata[pcadata >= 0], na.rm = T)
  pcadata[is.na(pcadata)] <- minValue/2
}

if (narm)
{
  minValue <- min(pcadata[pcadata >= 0], na.rm = T)
  pcadata[is.na(pcadata)] <- minValue/2
}


# More optional parameters:
# Log transform performed with the Generalized logarithm, like in metaboanalyst
# This method deals with zeros and negative numbers.
# The formula for base2 generalized log transform is: glog(x) = log2 ( x+sqrt(x^2 + a^2) ) / 2
# With the constant a = 1... this seems to be not true..
# Generalized Log transformation:
if (log == TRUE)
{
# tpcadata <- pcadata
# pcadata[tpcadata <= 0] <- log2( pcadata[tpcadata <= 0]/2 + sqrt(pcadata[tpcadata <= 0]^2+1)/2)
# pcadata[tpcadata > 0] <- log2(pcadata[tpcadata > 0])

 # pcadata <- sapply(pcadata, function(x) log2( x/2 + sqrt(x^2 + 1^2)/2   ))
  pcadata <- log2(pcadata)
# pcadata <- log2( pcadata/2 + sqrt(pcadata^2+.00000005)/2)

}

if (power == TRUE)
{
  pcadata <- sqrt(pcadata)
}

# Do the PCA, can choose between mean center, no scaling ,auto and pareto scaling.
# Need to add range, vast and level scaling,
# see: "Centering, scaling, and transformations: improving the biological information content of metabolomics data"

if (scale == "auto") {
  #pcadata <- scale(pcadata)
  pcadata <- scale(pcadata)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)
}

else if(scale == "pareto") {
  pcadata <- pareto(pcadata)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)
}

else if(scale == "none"){
  pcadata <- scale(pcadata, center = center, scale = F)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)
}

else if(scale == "range") {
  pcadata <- Range(pcadata)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)
}
else if(scale == "vast") {
  pcadata <- Vast(pcadata)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)
}
else if(scale == "level") {
  pcadata <- LevelScaling(pcadata)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)
}
else {
  cat("Incorrect type of scaline, data has not been scaled. \n")
  pcadata <- scale(pcadata, center = center, scale = F)
  met.pca <- prcomp(na.omit(pcadata),center = F, scale. = F)

}

# Here select which PCs to use, default 1 and 2
# naming system here is that PC1 refers to the x axis and PC2 to the y axis.
# Can be whichever PC selected.

# get data for percentage of variance per PC
# Gives the possibility of creating a scree plot.
pca.sum <- summary(met.pca)

PC1p <- pca.sum$importance[2,pcs[1]] * 100
PC2p <- pca.sum$importance[2,pcs[2]] * 100

PC1lab <- paste(c("PC "), pcs[1],c(" ("), round(PC1p,1), c("%)"), sep ="")
PC2lab <- paste(c("PC "), pcs[2],c(" ("), round(PC2p,1), c("%)"), sep ="")

# For the first iteration, simple plot:
Loadings <- as.data.frame(met.pca$rotation)


# calculate PC1:
nsamples <- length(pcadata[,1])
PC1 <- matrix(0, nsamples)
PC2 <- matrix(0, nsamples)

for (i in 1:nsamples)
  {
    PC1[i] <- sum(pcadata[i,]*Loadings[,pcs[1]])
    PC2[i] <- sum(pcadata[i,]*Loadings[,pcs[2]])
}



ggdata <- data.frame(catdata, PC1, PC2)
if (groups != "none")
{
  gro <- which(colnames(ggdata) == groups)
  ggdata[,gro] <- as.factor(ggdata[,gro])
  glabels <- NULL
  for (i in 1:nlevels(ggdata[,gro]))
    {
    b <- levels(ggdata[,gro])[i]
    glabels <- c(glabels,b)
    }
}



dashcol = c("snow4")

maxPC1 <- ceiling(max(ggdata$PC1))+4
maxPC2 <- ceiling(max(ggdata$PC2))+4
minPC1 <- floor(min(ggdata$PC1))-4
minPC2 <- floor(min(ggdata$PC2))-4

plotColors <- colors

p <- ggplot(data = ggdata, aes(x = PC1, PC2))+
      geom_hline(yintercept = 0, color = dashcol ,linetype = "dashed") +
      geom_vline(xintercept = 0, color = dashcol ,linetype = "dashed")

# Add the confidence interval
if (CI == TRUE) {
  p <- p + stat_ellipse(geom = "polygon",alpha = .3,
       aes(fill=ggdata[,gro], color = NA), show.legend = F, type = "norm", linetype = 3)
}


if (groups == "none")
  {
p <- p +  geom_point(aes(fill = plotColors[1]),shape = 21, size = pointsize)+
          scale_x_continuous(name = PC1lab) +
          scale_y_continuous(name = PC2lab)+
          theme_bw() +
          scale_fill_manual(values = plotColors, name = "", labels = "") +
          theme(legend.position = "none")
  }
else {
p <- p +  geom_point(aes(fill = ggdata[,gro]),shape = 21, size = pointsize)+
          scale_x_continuous(name = PC1lab) +
          scale_y_continuous(name = PC2lab)+
          scale_fill_manual(values=plotColors,name = "") +
          scale_color_manual(values=plotColors,name = "", labels = glabels) +

          theme_bw()
}



if (bLabel == TRUE)
  {
  Labels <- which(colnames(ggdata) == label)
  p <- p + geom_text_repel(aes(label=ggdata[,Labels]),
                           segment.color = 'black',
                           segment.size = 0.5,
                           arrow = arrow(length = unit(0.01, 'npc')),
                           force = 1)
  }


p <- p + ggtitle("Scores plot")+
  theme( plot.title = element_text(hjust = 0.5))


if (corrList == TRUE)
{

  corloadings1 <- cor(pcadata, met.pca$x)
  corvectors <- matrix(0,length(corloadings1[,1]), 1)

  metnames <- colnames(pcadata)
  corloadings <- data.frame(corloadings1,corvectors,metnames)
  mostpc <- sort(abs(corloadings1[,pcs[1]]), decreasing = TRUE)

  usethis <- subset(corloadings1, abs(corloadings1[,pcs[1]]) > mostpc[nh+1])

  print(row.names(usethis))
  outp = row.names(usethis)
}





if (corplot == TRUE)
  {

  corloadings1 <- cor(pcadata, met.pca$x)
  corvectors <- matrix(0,length(corloadings1[,1]), 1)

  for (i in 1:length(corvectors))
  {
  corvectors[i]  <- sqrt( ((corloadings1[i,pcs[1]])^2) + ((corloadings1[i,pcs[2]])^2))
  }

  metnames <- colnames(pcadata)
  corloadings <- data.frame(corloadings1,corvectors,metnames)
  most <- sort(abs(corvectors), decreasing = TRUE)


  usethis <- subset(corloadings, corloadings$corvectors > most[nh+1])

  a <- seq(0, 2*pi, length = 100)
  plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = PC1lab,  ylab = PC2lab)
  abline(h = 0, v = 0, lty = 2)
  points(corloadings1[,pcs[2]]~corloadings1[,pcs[1]])
  points(usethis[,pcs[2]]~ usethis[,pcs[1]], pch = 21,bg = "red", col = "red" )
  text(usethis[,pcs[2]]~usethis[,pcs[1]],labels = usethis[1:nh,ncol(usethis)],
     cex = 0.7, pos = 4, offset = 0.5)

}






if (scree == TRUE) {
  # draw the scree plot... accumulated and explaned...first 5 components.
  cumul <- rep(0,5)
  pccon <- rep(0,5)
  cumulL <- rep("",5)
  pcconL <- rep("",5)
  for (i in 1:5)
  {
    pccon[i] <- pca.sum$importance[2,i]
    cumul[i] <- sum(pccon)
    pcconL[i] <- paste(round((pccon[i]*100),1), "%", sep = "")
    cumulL[i] <- paste(round((cumul[i]*100),1), "%", sep = "")
  }

  plot(cumul~c(1:5), cex = 1.1,xlab = "PC index",ylab = "Variance explained",
       xlim = c(1,5.5), col = "red" , ylim = c(0,max(cumul)+0.15*max(cumul)), frame.plot = F)
  lines(cumul~c(1:5), col = "green")
  points(pccon~c(1:5), cex = 1.1, col = "red")
  lines(pccon~c(1:5), col = "blue")
  abline(v = c(1:5))
  text(cumul~c(1:5), labels = cumulL, adj = c(-.2,0),srt = 45)
  text(pccon~c(1:5), labels = pcconL, adj = c(-.2,0),srt = 45)
}





if (loadinglist == TRUE)
{
  output <- as.data.frame(met.pca$rotation)
  row.names(output) <- t(colnames(pcadata))
}
else if (corrList == TRUE) {
  return(outp)
}
else {
  output <- p
}

return(output)

}





