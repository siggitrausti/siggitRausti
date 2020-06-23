#' # A quite raw function which can find ratios that are more interesting than single pairs of predictors in lin.reg.models
#'
#' @param dataset,covariates,start,p_thresh
#' @keywords data_processing
#' @export
#' @examples
#' p_gain_analysis()

p_gain_analysis <- function(dataset,covariates,start,p_thresh){
  # The structure of this dataset assumes that the first column is the predicted variable and that this is (for now)
  # just applicable for logistic regression modeling......
  
  
  require(gridExtra)
  require(siggitRausti)
  require(ggplot2)
  #require(arm)

  if(missing(p_thresh)){
    p_thresh <- 0.05
  }
  houston_fin_modules <- dataset
  colnames(houston_fin_modules)[1] <- 'Death'
  for (i in start:ncol(houston_fin_modules)){
    houston_fin_modules[,i] <- min_max(houston_fin_modules[,i],0.1,1)
  }
  
  if (is.factor(houston_fin_modules[,1])){
    p_vec <- rep(NA,ncol(houston_fin_modules))
    for (i in start:ncol(houston_fin_modules)){
      if (start != 2){
        covariates_temp <- c(covariates,colnames(houston_fin_modules)[i])
        formula_if_covariates <- reformulate(covariates_temp,'Death')
        temp_model <- glm(formula_if_covariates,data = houston_fin_modules,family = binomial(link = "logit"))
      } else {
        temp_model <- glm(Death ~ houston_fin_modules[,i],data = houston_fin_modules,family = binomial(link = "logit"))
      }
      p_vec[i] <- coef(summary(temp_model))[,4][start]
    }
    #p_vec <- p.adjust(p_vec,method = 'BH')
    id_sign <- which(p_vec < 0.05)
    if (length(id_sign) == 0){
      stop("ERROR: No feature is significantly correlated with the trait")
    }
    id_lowest <- which(p_vec == min(p_vec,na.rm=T))
    
    
    # Create a ratio vector:
    # First, create a dataframe with all the ratios:
    res <- do.call(cbind, combn(houston_fin_modules[,c(start:ncol(houston_fin_modules))], 2, FUN= function(x) list(x[2]/x[1])))
    names(res) <- combn(names(houston_fin_modules[,c(start:ncol(houston_fin_modules))]), 2, FUN = function(x) paste(x[2], x[1], sep="_"))
    res <- res[,is.finite(colSums(res))]
    
    houston_fin_modules_1 <- data.frame(cbind(houston_fin_modules[,c(1:start-1)],res))
    colnames(houston_fin_modules_1)[1] <- 'Death'
    #colnames(houston_fin_modules3)[ncol(houston_fin_modules3)] <- combn(names(houston_fin_modules_SUBSET[,-c(1,2)]), 2, FUN = function(x) paste(x[2], x[1], sep="_"))
    #houston_fin_modules3 <- houston_fin_modules3[,is.finite(colSums(houston_fin_modules3))]
    
    keep_vec <- c()
    value_to_contain <- colnames(houston_fin_modules)[id_sign]
    for (i in start:ncol(houston_fin_modules_1)){
      values <- strsplit(colnames(houston_fin_modules_1)[i],'_')[[1]]
      if (any(value_to_contain %in% values)){
        keep_vec <- c(keep_vec,i)
      }
    }
    
    houston_fin_modules_1 <- houston_fin_modules_1[,c(1:start-1,keep_vec)]
    
    # LOG TRANSFORM FOR REVERSIBILITY
    #houston_fin_modules_1[,(start:ncol(houston_fin_modules_1))] <- log(houston_fin_modules_1[,(start:ncol(houston_fin_modules_1))],2)
    
    
    p_vec_ratios <- rep(NA,ncol(houston_fin_modules_1))
    for (i in start:ncol(houston_fin_modules_1)){
      if (start != 2){
        covariates_temp <- c(covariates,colnames(houston_fin_modules_1)[i])
        formula_if_covariates <- reformulate(covariates_temp,'Death')
        #print(formula_if_covariates)
        temp_model <- glm(formula_if_covariates,data = houston_fin_modules_1,family = binomial(link = "logit"))
      } else {
        temp_model <- glm(Death ~ houston_fin_modules_1[,i],data = houston_fin_modules_1,family = binomial(link = "logit"))
      }
      p_vec_ratios[i] <- coef(summary(temp_model))[,4][start]
    }
    #p_vec_ratios <- p.adjust(p_vec_ratios,method = 'BH')
    
    
    # Now calculate the p-gain:
    p_gain <- rep(NA,ncol(houston_fin_modules_1))
    for (i in start:ncol(houston_fin_modules_1)){
      values <- strsplit(colnames(houston_fin_modules_1)[i],'_')[[1]]
      id_1 <- which(colnames(houston_fin_modules) %in% values[1])
      p_1 <- p_vec[id_1]
      id_2 <- which(colnames(houston_fin_modules) %in% values[2])
      p_2 <- p_vec[id_2]
      if (p_2 < p_1){
        p_gain[i] = p_2/p_vec_ratios[i]
      } else if (p_2 > p_1){
        p_gain[i] = p_1/p_vec_ratios[i]
      }
    }
    #plot(sort(p_gain))
    
    #id_sign <- which(p_gain > 1/(2*p_thresh)) # No correction for multiple comparisons....
    #print(max(p_gain,na.rm=T))
    id_sign2 <- which(p_gain == max(p_gain,na.rm=T)) # No correction for multiple comparisons....
    #print(id_sign)
    
    if (p_gain[id_sign2] > (1/(2*p_thresh))){
      print(paste0('Significant improvement in ratios in ',colnames(houston_fin_modules_1)[id_sign2]))
      houston_outp <- houston_fin_modules_1[,c(1:start-1,id_sign2)]
      colnames(houston_outp)[1] <- colnames(dataset)[1]
      print(paste0('Best single gene predictor: ',colnames(houston_fin_modules)[id_lowest]))
      houston_outp_single <- houston_fin_modules[,c(1:start-1,id_lowest)]
      #return(houston_outp)
      source('./ROC_maker.R') # Clearly need to turn this into a function within siggitRausti
      ROC_single = ROC_maker(houston_outp_single[complete.cases(houston_outp_single),]) + ggtitle(paste0(colnames(houston_fin_modules)[id_lowest]))
      ROC_ratio = ROC_maker(houston_outp[complete.cases(houston_outp),]) + ggtitle(paste0('Ratio of ',colnames(houston_fin_modules_1)[id_sign2]))
      #p_gain <- p_gain[!is.na(p_gain)]
      data_pgain <- data.frame(cbind(Sequence = seq(1,length(p_gain),by = 1),Gain = p_gain))
      data_pgain$NAME <- colnames(houston_fin_modules_1)
      data_pgain <-data_pgain[order(data_pgain$Gain,decreasing = F),]
      data_pgain$Sequence2 <- seq(1,nrow(data_pgain),by=1)
      #data_pgain <- data_pgain[complete.cases(data_pgain),]
      #print(data_pgain)
      pgain_image = ggplot(data_pgain,aes(x=Sequence2,y=Gain,label = NAME)) + geom_point() + 
        geom_text(aes(label=ifelse(Gain>(1/(2*p_thresh)),as.character(NAME),'')),size = 6,position = position_nudge(y = -0.05*max(p_gain,na.rm=T),x=-0.2*length(p_gain)))
      pgain_image <- plotLookForPaper(pgain_image,'Gain value','') + geom_hline(yintercept=1/(2*p_thresh), linetype="dashed", color = "red") + 
        xlim(0,max(data_pgain$Sequence2,na.rm=T)+1)
      grid.arrange(pgain_image, ROC_single, ROC_ratio, ncol=2, nrow =2)
    } else {
      print('No significant improvements in ratios over single metabolites')
      if (length(id_sign) != 0){
        houston_outp <- houston_fin_modules[,c(1:start-1,id_lowest)]
        colnames(houston_outp)[1] <- colnames(dataset)[1]
        single_image <- ROC_maker(houston_outp[complete.cases(houston_outp),]) + ggtitle(paste0(colnames(houston_fin_modules)[id_lowest]))
        grid.arrange(single_image, ncol=1, nrow =1)
      } else {
        houston_outp <- NULL
      }
    }
  } else {
    p_vec <- rep(NA,ncol(houston_fin_modules))
    for (i in start:ncol(houston_fin_modules)){
      if (start != 2){
        covariates_temp <- c(covariates,colnames(houston_fin_modules)[i])
        formula_if_covariates <- reformulate(covariates_temp,'Death')
        temp_model <- lm(formula_if_covariates,data = houston_fin_modules)
      } else {
        temp_model <- lm(Death ~ houston_fin_modules[,i],data = houston_fin_modules)
      }
      p_vec[i] <- coef(summary(temp_model))[,4][start]
    }
    #p_vec <- p.adjust(p_vec,method = 'BH')
    id_sign <- which(p_vec < 0.05)
    id_lowest <- which(p_vec == min(p_vec,na.rm=T))
    
    
    # Create a ratio vector:
    # First, create a dataframe with all the ratios:
    res <- do.call(cbind, combn(houston_fin_modules[,c(start:ncol(houston_fin_modules))], 2, FUN= function(x) list(x[2]/x[1])))
    names(res) <- combn(names(houston_fin_modules[,c(start:ncol(houston_fin_modules))]), 2, FUN = function(x) paste(x[2], x[1], sep="_"))
    res <- res[,is.finite(colSums(res))]
    
    houston_fin_modules_1 <- data.frame(cbind(houston_fin_modules[,c(1:start-1)],res))
    colnames(houston_fin_modules_1)[1] <- 'Death'
    #colnames(houston_fin_modules3)[ncol(houston_fin_modules3)] <- combn(names(houston_fin_modules_SUBSET[,-c(1,2)]), 2, FUN = function(x) paste(x[2], x[1], sep="_"))
    #houston_fin_modules3 <- houston_fin_modules3[,is.finite(colSums(houston_fin_modules3))]
    
    keep_vec <- c()
    value_to_contain <- colnames(houston_fin_modules)[id_sign]
    for (i in start:ncol(houston_fin_modules_1)){
      values <- strsplit(colnames(houston_fin_modules_1)[i],'_')[[1]]
      if (any(value_to_contain %in% values)){
        keep_vec <- c(keep_vec,i)
      }
    }
    
    houston_fin_modules_1 <- houston_fin_modules_1[,c(1:start-1,keep_vec)]
    
    # LOG TRANSFORM FOR REVERSIBILITY
    #houston_fin_modules_1[,(start:ncol(houston_fin_modules_1))] <- log(houston_fin_modules_1[,(start:ncol(houston_fin_modules_1))],2)
    
    p_vec_ratios <- rep(NA,ncol(houston_fin_modules_1))
    for (i in start:ncol(houston_fin_modules_1)){
      if (start != 2){
        covariates_temp <- c(covariates,colnames(houston_fin_modules_1)[i])
        formula_if_covariates <- reformulate(covariates_temp,'Death')
        temp_model <- lm(formula_if_covariates,data = houston_fin_modules_1)
      } else {
        temp_model <- lm(Death ~ houston_fin_modules_1[,i],data = houston_fin_modules_1)
      }
      p_vec_ratios[i] <- coef(summary(temp_model))[,4][start]
    }
    #p_vec_ratios <- p.adjust(p_vec_ratios,method = 'BH')
    #print(p_vec_ratios)
    
    # Now calculate the p-gain:
    p_gain <- rep(NA,ncol(houston_fin_modules_1))
    for (i in start:ncol(houston_fin_modules_1)){
      values <- strsplit(colnames(houston_fin_modules_1)[i],'_')[[1]]
      id_1 <- which(colnames(houston_fin_modules) %in% values[1])
      p_1 <- p_vec[id_1]
      id_2 <- which(colnames(houston_fin_modules) %in% values[2])
      p_2 <- p_vec[id_2]
      if (p_2 < p_1){
        p_gain[i] = p_2/p_vec_ratios[i]
      } else if (p_2 > p_1){
        p_gain[i] = p_1/p_vec_ratios[i]
      }
    }
    #plot(sort(p_gain))
    
    #id_sign <- which(p_gain > 1/(2*p_thresh)) # No correction for multiple comparisons....
    #print(max(p_gain,na.rm=T))
    id_sign2 <- which(p_gain == max(p_gain,na.rm=T)) # No correction for multiple comparisons....
    #print(id_sign)
    
    if (p_gain[id_sign2] > (1/(2*p_thresh))){
      print(paste0('Significant improvement in ratios in ',colnames(houston_fin_modules_1)[id_sign2]))
      houston_outp <- houston_fin_modules_1[,c(1:start-1,id_sign2)]
      colnames(houston_outp)[1] <- colnames(dataset)[1]
      print(paste0('Best single gene predictor: ',colnames(houston_fin_modules)[id_lowest]))
      houston_outp_single <- houston_fin_modules[,c(1:start-1,id_lowest)]
      #return(houston_outp)
      #source('./ROC_maker.R') # Clearly need to turn this into a function within siggitRausti
      #ROC_single = ROC_maker(houston_outp_single[complete.cases(houston_outp_single),]) + ggtitle(paste0(colnames(houston_fin_modules)[id_lowest]))
      #ROC_ratio = ROC_maker(houston_outp[complete.cases(houston_outp),]) + ggtitle(paste0('Ratio of ',colnames(houston_fin_modules_1)[id_sign2]))
      #p_gain <- p_gain[!is.na(p_gain)]
      data_pgain <- data.frame(cbind(Sequence = seq(1,length(p_gain),by = 1),Gain = p_gain))
      data_pgain$NAME <- colnames(houston_fin_modules_1)
      data_pgain <-data_pgain[order(data_pgain$Gain,decreasing = F),]
      data_pgain$Sequence2 <- seq(1,nrow(data_pgain),by=1)
      #data_pgain <- data_pgain[complete.cases(data_pgain),]
      #print(data_pgain)
      pgain_image = ggplot(data_pgain,aes(x=Sequence2,y=Gain,label = NAME)) + geom_point() + 
        geom_text(aes(label=ifelse(Gain>(1/(2*p_thresh)),as.character(NAME),'')),size = 6,position = position_nudge(y = -0.05*max(p_gain,na.rm=T),x=-0.2*length(p_gain)))
      pgain_image <- plotLookForPaper(pgain_image,'Gain value','') + geom_hline(yintercept=1/(2*p_thresh), linetype="dashed", color = "red") + 
        xlim(0,max(data_pgain$Sequence2,na.rm=T)+1)
      grid.arrange(pgain_image, ncol=1, nrow =1)
    } else {
      print('No significant improvements in ratios over single metabolites')
      if (length(id_sign) != 0){
        houston_outp <- houston_fin_modules[,c(1:start-1,id_lowest)]
        colnames(houston_outp)[1] <- colnames(dataset)[1]
        #single_image <- ROC_maker(houston_outp[complete.cases(houston_outp),]) + ggtitle(paste0(colnames(houston_fin_modules)[id_lowest]))
        #grid.arrange(single_image, ncol=1, nrow =1)
      } else {
        houston_outp <- NULL
      }
    }
    covariates_single <- c(covariates,colnames(houston_fin_modules)[id_lowest])
    lm_single <- lm(reformulate(covariates_single,'Death'),data = houston_fin_modules)
    covariates_ratio <- c(covariates,colnames(houston_fin_modules_1)[id_sign2])
    lm_ratio <- lm(reformulate(covariates_ratio,'Death'),data = houston_fin_modules_1)
    print(summary(lm_single))
    print(summary(lm_ratio))
    
  }
  return(houston_outp)
}