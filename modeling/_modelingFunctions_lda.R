# =====================================================#
# @author: Arnau Serra-Cayuela
# @contact: serracay@ualberta.ca
# @group: WishartLab (University of Alberta)
# @date: 2017
# 

# ====================================================#

### --- Set the working directory
#setwd("/media/data/statAnalysis/MYCO/GWAS/")

### --- Load libraries

library(usdm) # vif to check multicollinearity
library(mlr)
library(jsonlite)
source("modeling/_splitIndependentDataSet.R")
source("modeling/_summaryFunctions.R")
require(pROC)


### ==  functions ==

### --- Models (classifier and maybe feature selection)

## -- Logistic regression
ols <- function(task, outerFolds, outerRep, seed){
  # Generate the learner
  classifier <- "classif.logreg"
  lrn <- mlr::makeLearner(classifier, predict.type = "prob")
  ## Create cross-validation schema
  outer <- mlr::makeResampleDesc("RepCV", folds = outerFolds, reps = outerRep )
  # ## estimate performance in the train subset
  set.seed(seed)
  r <- mlr::resample(lrn, task, resampling = outer, show.info = TRUE, measures = list(mlr::auc, mlr::tpr, mlr::tnr, mlr::acc))
  return(list(r, lrn))
  #return(r)
}


## -- LDA
lda <- function(task, outerFolds, outerRep, seed){
  # Generate the learner
  classifier <- "classif.lda"
  lrn <- mlr::makeLearner(classifier, predict.type = "prob")
  ## Create cross-validation schema
  outer <- mlr::makeResampleDesc("RepCV", folds = outerFolds, reps = outerRep )
  # ## estimate performance in the train subset
  set.seed(seed)
  r <- mlr::resample(lrn, task, resampling = outer, show.info = TRUE, measures = list(mlr::auc, mlr::tpr, mlr::tnr, mlr::acc))
  return(list(r, lrn))
  #return(r)
}



## -- Ridge regression (penalized logistic regression)
ridgeReg <- function(task, outerFolds, outerRep, seed){
  # Generate the learner
  classifier <- "classif.glmnet"
  lrn <- mlr::makeLearner(classifier, predict.type = "prob")
  # define parameters we want to tune 
  ps <- makeParamSet( makeDiscreteParam("s", values = seq(0, 1, by=0.1)), # reduced from by = 0.1, to be faster
                      makeDiscreteParam("alpha", values = 0))
  
  # slowest method, but it gives all the possibilities
  ctrl <- makeTuneControlGrid()
  # alternatives
  #ctrl <- makeTuneControlRandom(maxit = 6L, same.resampling.instance = FALSE)
  #ctrl <- makeTuneControlIrace(maxExperiments = 200L)
  #ctrl <- makeTuneControlGenSA(budget = 100L) 
  
  # do a nested 5 fold cross-validation
  inner <- makeResampleDesc("CV", iters = 3L)
  # wrap the tunable parameters in the learner
  lrn <- mlr::makeTuneWrapper(lrn, resampling = inner, par.set = ps,
                              control = ctrl, show.info = TRUE, measures = list(mlr::auc))
  
  ## Create outer cross-validation schema
  outer <- mlr::makeResampleDesc("RepCV", folds = outerFolds, reps = outerRep )
  ## estimate performance in the train subset
  set.seed(seed)
  ## do CV
  r <- mlr::resample(lrn, task, resampling = outer, show.info = TRUE, extract = getTuneResult,
                     measures = list(mlr::auc, mlr::tpr, mlr::tnr, mlr::acc))
  return(list(r, lrn))
}


### --- personalized ROC plot function
makeROCplot <- function(predicted, legend = TRUE){
  ro <- pROC::roc(predicted$data$truth, predicted$data$prob.control) # before was prob.case
  par(mfrow=c(1,1))
  sen = ro$sensitivities
  spe = 1 - ro$specificities
  plot(sen~spe, type = "n", las=1, xlab = "1-Specificity (False positive rate)",
       ylab = "Sensitivity (True positive rate)")
  abline(a = 0, b = 1, col = "grey80", lty=2)
  #grid(col = "lightgray", lty = "dotted", lwd = 1)
  points(sen~spe, type = "l", lw=2, col = "#4911e7")
  if(legend){
    cint <- ci(ro)
    text(0.7, 0.2, sprintf("AUC (95%% CI) = %1.2f (%1.2f-%1.2f)" , cint[2], cint[1], cint[3]))
  }
  saveRDS(ro, "ro.rds")
  print("ROC analysis done!")
}

#example
# makeROCplot(pred)
# ro <- pROC::roc(pred$data$truth, pred$data$prob.case)
# mlr::performance(pred, measures = list(mlr::auc))



### --- Save model coefficient to JSON file
saveCoefJSON <- function(coef, classifier, exitPath){
  if(classifier == "classif.glmnet"){
    coefDF <- data.frame(SNP = names(coef[,1]), coef = coef[,1])
  } else{
    coefDF <- data.frame(SNP = names(coef), coef = coef)
  }
  rownames(coefDF) <- NULL
  coefJSON <- toJSON(coefDF)#, pretty = TRUE)
  write(coefJSON, paste(exitPath, "coef.json", sep="/"))
}


## Main function, which estimate the coefficients

doAnalysis <- function(datasetPath, exitPath, groupingVar = "class",  positiveResponse = "case",
                       pSplit = 0.8,  #classifier = "classif.logreg", 
                       outerFolds = 2L,
                       outerRep = 2L, seed = 1984, verbose = FALSE){
  
  ## remove previous warnings
  #assign("last.warning", NULL, envir = baseenv())
  
  ## load data
  datas <- read.csv(datasetPath)
  da <- data.frame(datas[,-1], row.names=datas[,1])
  ## study dimensions
  if(verbose){
    cat("\nData Set dimensions: ", dim(da))
  }
  
  info <- list()
  info$nSNP <- dim(da)[2] - 1 # first column is the grouping 
  info$indv <- summary(da[,1])
  
  
  ## get the study name
  study <- gsub(".*\\/","", datasetPath)
  study <- gsub("_Resampled.csv", "", study)
  ## get the study name/number
  studyID <- sub('.*\\_', '', study)
  info$studyID <- studyID
  ## name of the disease
  diseaseName <- sub("_+[^_]*$", "", study)
  #diseaseName <- gsub("_"," ", diseaseName)
  diseaseName <- gsub("_+", " ", diseaseName)
  diseaseName <- gsub("  ", " ", diseaseName)
  
  capitalize <- function (string) 
  {
    capped <- grep("^[A-Z]", string, invert = TRUE)
    substr(string[capped], 1, 1) <- toupper(substr(string[capped], 
                                                   1, 1))
    return(string)
  }
  info$diseaseName  <- capitalize(diseaseName)
 
  
  ## optional, name of the dataset
  # if(exists("studyName")){
  #   datasetName <- studyName
  # } else{
  #   datasetName <- "GWAS"
  # }
  
  
  ## stratifying variable
  stratifyBy <- groupingVar
  
  ## == split data into training and testing
  dataSplitIndep_idx <- splitData_idx(da, strat = stratifyBy, pSplit = pSplit,
                                      seed = seed)
  ## create sub datasets
  train.set <- da[dataSplitIndep_idx[[1]],]
  test.set <- da[-dataSplitIndep_idx[[1]],]
  if(verbose){
    cat("\nData split in train and test\n")
  }
  
  
  ## === check for multicollinarity and assign method
  mcoli <- usdm::vif(datas[,-1])
  if(verbose){
    cat("\nEstimated VIF\n")
    print(mcoli)
  }
  
  ## Define the task
  task <- mlr::makeClassifTask(data = train.set, target = groupingVar, positive = positiveResponse)
  
  ## choose the model to use depending on the mcoli output, and resample
  # if(classifier2 == "LDA"){
  #   cat("\nldaaaaaaaaaaaaaa\n")
  #   method <- "lda"
  #   classifier <- "classif.lda"
  #   r<-lda(task, outerFolds, outerRep, seed)
  # } else{
  
    # if(classifier == "classif.glmnet"){
    #   cat("\nPenalized regression\n")
    #   method <- "penalized"
    #   classifier <- "classif.glmnet"
    #   r <- ridgeReg(task, outerFolds, outerRep, seed)
    # } else if(classifier == "classif.logreg"){
    #     cat("\nOrdinary Least Squares\n")
    #     method <- "OLS"
    #     classifier <- "classif.logreg"
    #     r<-ols(task, outerFolds, outerRep, seed)
    # } else if(classifier == "classif.lda"){
    #     cat("\nldaaaaaaaaaaaaaa\n")
    #     
    #     method <- "lda"
    #     classifier <- "classif.lda"
    #     r <- lda(task, outerFolds, outerRep, seed)
    # }
    # 
    # info$method <- method
  
  if(sum(is.infinite(mcoli[,2][2]))!=0){
    cat("\nPenalized regression\n")
    method <- "penalized"
    classifier <- "classif.glmnet"
    r <- ridgeReg(task, outerFolds, outerRep, seed)
  } else{
    cat("\nOrdinary Least Squares\n")
    method <- "OLS"
    classifier <- "classif.logreg"
    r<-ols(task, outerFolds, outerRep, seed)
  }
  
  info$method <- method
  
  ## extract the values
  lrn <- r[[2]]
  r <- r[[1]]
  
  if(verbose){
    cat("\nresampling DONE!!!\n")
    print(r$measures.test)
  }
  
  ## == plot metriccs results
  png(filename = paste(exitPath, "trainingPlotMetrics.png", sep="/"), width = 2000, height = 2000, res = 260)
  plotTrainingMetrics(r)
  dev.off()
  ## show as table
  trainingMetrics <- tableMetricTrain(r, asDF = TRUE); trainingMetrics
  # save as a csv
  print("Writing training data")
  write.csv(trainingMetrics, paste(exitPath, "trainingPerformance.csv", sep="/"))
  
  
  ### === Build the final model
  
  if(method == 'lda'){
    print("lda")
    modWarn <- tryCatch.W.E(
      mod <- mlr::train(lrn, task)
    )
  }
  if(method == "OLS"){
    print("OLS")
    modWarn <- tryCatch.W.E(
      mod <- mlr::train(lrn, task)
    )
    coef <- coef(mod$learner.model)
    print(coef)
    # save coef as a JSON file
    # saveCoefJSON(coef, classifier, exitPath) 
    write.csv(coef, paste(exitPath, "coef.csv", sep="/"), col.names = F)
    
    
  } else if(method =="penalized"){
    print("Building penalized model")
    # get the best tuned values and build the model
    tuneRes_S <- getNestedTuneResultsX(r)$s
    tuneRes_S_idx <- which(table(tuneRes_S) == max(table(tuneRes_S)))
    best_S <-tuneRes_S[tuneRes_S_idx][1] ## [1] in case there is more than one possible best S
    if(verbose){
      cat("\nBest S:\n" , best_S)
    }
    
    ## add the tuned parameters
    lrn <- mlr::makeLearner(classifier, predict.type = "prob")
    lrn <- setHyperPars(lrn, par.vals = list(s=best_S, alpha = 0))
    
    ## build the model
    modWarn <- tryCatch(
      mod <- mlr::train(lrn, task),
      warning=function(w) print(paste(w,"Model may be wrong, perfect fit", sep=":::")))
    
    # save coef as a JSON file
    #coef <- coef(mod$learner.model)
    #saveCoefJSON(coef, classifier, exitPath) 
    ### TODO: export the best parameters
    
    
  }
  
  info$warn <- modWarn
  
  ## === Predict with testing subset
  
  ## -- prediction of test subset
  print(mod)
  pred <- predict(mod, newdata = test.set)
  ## -- save into a table
  testMetrics <- performance(pred, measures = list(mlr::auc, mlr::tpr, mlr::tnr, mlr::acc))
  testMetrics <- as.data.frame(testMetrics)
  rownames(testMetrics) <- c("AUC", "Sensitivity", "Specificity", "Accuracy")
  print(testMetrics)
  write.csv(testMetrics, paste(exitPath, "testPerformance.csv", sep="/"))
  
  ## plot ROC curve
  png(filename = paste(exitPath, "ROCplot.png", sep="/"), width = 2000, height = 2000, res = 260)
  makeROCplot(pred)
  dev.off()

  ## return values
  trainRes <- list(outerCV = r$measures.test, summary = trainingMetrics)
  allRes <<- list(training = trainRes, model = mod, prediction = pred, test = testMetrics, info = info)
  return(list(allRes=allRes, method=method))
}


## examples
# dPath <- "data/resampled/out/Age_related_macular_degeneration_HGVST1565_Resampled.csv"
# analRes <- doAnalysis(dPath, verbose = TRUE)
# analRes


### Function to get potential warnings
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


### --- function to analyze all the studie in a specified folder
batchAnalysis <- function(inPath, outPath, studies, ...){
  ## -- function to change the name 
  pcName <- function(nom){
    return(gsub('[-? (),\'%]', '_', nom))
  }
  
  #studies <- dir(inPath)
  ## only csv files
  #studies <- studies[grep(".csv", studies)]
  
  #studies = c("Combining Information from Common Type 2 Diabetes Risk Polymorphisms Improves Disease Prediction_Resampled.csv")
  
  t1 <- Sys.time()
  
  for(i in studies){
    ## -- create a folder for each study
    fname <- pcName(i)
    fname <- substr(fname,1, nchar(fname)-14) # remove _resampled..
    print("-------------------------------------------------------------------")
    print(paste("Analizing study:", fname, sep=" "))
    print("-------------------------------------------------------------------")
    ## create a folder for this study
    dir.create(paste(outPath, fname, sep="/"))
    ## -- set the path to the study
    dPath <- paste(inPath, i, sep="/")
    ## -- output results
    exitPath <- paste(outPath, fname, sep="/")
    ## do the analysis
    analRes <- doAnalysis(dPath, exitPath, ...)
    ## add study info
    analRes$info$studyName <- i
    ## save R object
    saveRDS(analRes, file = paste(outPath, fname, "analysis.rds", sep="/"))
    saveRDS(analRes$info, file = paste(outPath, fname, "info.rds", sep="/"))
    aucanalysis <<- analRes$allRes[[4]][[1,1]]
    regressiontype <- analRes$method
    print("AUC: ")
    print(aucanalysis)
    return(list(aucanalysis=aucanalysis, regressiontype=regressiontype))

  }
  print("DONE! :)")
  print(paste("Total time required:", sprintf("%2.2f",Sys.time()-t1), "seconds", sep=" "))
  
}
