# =====================================================#
# @author: Arnau Serra-Cayuela
# @contact: serracay@ualberta.ca
# @group: WishartLab (University of Alberta)
# @date: 2017
# 

# ====================================================#

source("modeling/_stratification.R")


# dataset: 
# strat: Indicate the name of the strtifying variable (if any)
# pSplit : Percentage [0,1] of observations that may have the modeling subset
# seed: number for reproducible values

splitData_idx <- function(dataset, strat = FALSE, pSplit, seed){
  # number of observations
  n = dim(dataset)[1]
  if(strat!=FALSE){
    # add reference vector
    dataset$obsRef <- 1:n
    set.seed(seed)
    train.set_idx <- stratified(dataset, strat, pSplit)$obsRef
    test.set_idx <- setdiff(1:n, train.set_idx)
    return(list(train.set_idx, test.set_idx))
  } else{
    # in case non startification is required
    train.set_idx <- sample(n, size = pSplit*n)
    test.set_idx <- setdiff(1:n, train.set_idx)
    return(list(train.set_idx, test.set_idx))
  }
}




# pSplit = 0.8
# train.set = sample(n, size = 0.8*n)
# test.set = setdiff(1:n, train.set)

# library(caret)
# inTrain = createDataPartition(da$Species, p = 2/3, list = FALSE)
# dfTrain=df[inTrain,]
# dfTest=df[-inTrain,]
