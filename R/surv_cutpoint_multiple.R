#' # A function to check whether multiple genes can be used to discriminate groups in terms of survival curves
#' @param meta_data,features,time,event,minprop,format_signature
#' @keywords survival_analysis
#' @export
#' @examples
#' surv_cutpoint_multiple()


surv_cutpoint_multiple <- function(meta_data,features,time,event,minprop,format_signature){
  # A function to check whether multiple genes can be used to discriminate groups in terms of survival curves.
  # This function 1) scales data using pareto scaling and 2) z-scores data to come up with optimal signature cutoff. 
  paretoscale <- function(vector){
    x <- vector
    x.centered <- x - mean(x)
    # Then we perform scaling on the mean-centered matrix
    x.sc <- x.centered/sqrt(sd(x))
    return(x.sc)
  }
  
  
  require(survminer)
  if(missing(format_signature)){
    format_signature = rep('Harmful!',length(features))
  }
  z_list <- data.frame(matrix(data=NA,nrow=nrow(meta_data),ncol=length(features)))
  for (i in 1:length(features)){
    id_in_data <- which(colnames(meta_data) %in% features[i])
    pareto_scaled_dat <- paretoscale(meta_data[,id_in_data])
    z_scored <- as.numeric(scale(pareto_scaled_dat))
    z_list[,i] <- z_scored
  }
  #print(z_list)
  
  for (i in 1:length(format_signature)){
    if (format_signature[i] == 'Beneficial!'){
      z_list[,i] = -z_list[,i]
    }
  }
  #print(z_list)
  value_vec = rowSums(z_list)
  #print(value_vec)
  meta_data$value_vec = value_vec
  if (missing(minprop)){
    res.cut <- surv_cutpoint(meta_data, time = time, event = event,
                             variables = c('value_vec'))
  } else {
    res.cut <- surv_cutpoint(meta_data, time = time, event = event,
                             variables = c('value_vec'),minprop = minprop)
  }
  
  res.cat <- surv_categorize(res.cut)
  print(head(res.cat))
  fit <- survfit(Surv(survival, status) ~ value_vec, data = res.cat)
  p3 <- ggsurvplot(fit, data = res.cat, risk.table=F,pval=T, conf.int=F,
                   xlim = c(0,60),xlab = 'Time in months', # Choose length of follow-up time (120 = 10 years, 60 = 5 years)
                   ylab = 'Overall survival',
                   break.time.by = 12,
                   legend.labs = c('High signature','Low signature'),
                   #legend.labs = c('High ASL','Low ASL'),
                   palette = c('darkred','gray10'),
                   ggtheme = theme_bw(),
                   legend=c(0.3,0.15),
                   pval.size = 8,
                   pval.coord = c(3,0.375))
  p3 <- ggpar(p3,font.x = c(20, "bold", "gray15"),
              font.y = c(20, "bold", "gray15"),font.legend = c(14),font.tickslab = c(16)) 
  return(p3)
}