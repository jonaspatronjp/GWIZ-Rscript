# =====================================================#
# @author: Arnau Serra-Cayuela
# @contact: serracay@ualberta.ca
# @group: WishartLab (University of Alberta)
# @date: 2017
# 

# ====================================================#

createMetricText <- function(CVobject, metric, toDF = asDF){
  
  # estimate specificity
  # if(metric == "fpr"){
  #   CVobject$measures.test[,metric] <- 1 - CVobject$measures.test[,metric]
  # }
  
  if(toDF){
    # estimate CI
    n <- length(CVobject$measures.test[,metric])
    metricCI <- qt(0.975,df = (n-1))*sd(CVobject$measures.test[,metric], na.rm = TRUE)/sqrt(length(CVobject$measures.test[,metric]))
    metricMean <- mean(CVobject$measures.test[,metric], na.rm = TRUE)
    metricArray <- c( metricMean, metricMean - metricCI, metricMean + metricCI)
    return(metricArray)
    
  } else{
    # estimate CI
    n <- length(CVobject$measures.test[[metric]])
    metricCI<- qt(0.975,df = (n-1))*sd(CVobject$measures.test[[metric]], na.rm = TRUE)/sqrt(length(CVobject$measures.test[[metric]]))
    metricMean <- mean(CVobject$measures.test[[metric]], na.rm = TRUE)
    metricText <- paste( sprintf("%.3f", metricMean), " (", sprintf("%.3f", metricMean - metricCI), " - ", sprintf("%.3f", metricMean + metricCI), ")", sep = "" )
    return(metricText)
  }
} 


## format/style the output
tableMetricTrain <- function(CVobject, asDF = FALSE){
  
  if(asDF){
    metrics <- names(CVobject$measures.test)
    output <- matrix(ncol  = 3, nrow =length(metrics)-1)
    for(i in 2:length(metrics)){
      output[i-1,] <- createMetricText(CVobject, metrics[i], toDF = TRUE)
    }
    
    df <- as.data.frame(output)
    rownames(df) <- c("AUC", "Sensitivity", "Specificity", "Accuracy")
    names(df) <- c("mean", "low95%", "high95%")
    return(df)
  } else{
    df<- data.frame("AUC" = createMetricText(CVobject, "AUC"), "Sensitivity" = createMetricText(CVobject, "Sensitivity"),
                    "Specificity" = createMetricText(CVobject, "Specificity"), "Accuracy" = createMetricText(CVobject, "Accuracy"))
    return(df)
  }
}


plotTrainingMetrics <- function(CVobject){
  # estimate specificity
  #fpr_idx <- which(names(CVobject$measures.test) == "fpr")
  #CVobject$measures.test[[fpr_idx]] <- 1 - CVobject$measures.test[[fpr_idx]]
  
  names(CVobject$measures.test) <- c("iter", "AUC", "Sensitivity", "Specificity", "Accuracy")
  boxplot(CVobject$measures.test[,-1], las = 1, ylim = c(0,1), main = "Training values")
  abline(h = 0.5 , col = "grey80", lty = 2)
  
}

