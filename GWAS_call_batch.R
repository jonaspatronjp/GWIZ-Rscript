# Start of code

setwd("/Users/wishart/Desktop/GWAS_project_trimmed/")  ### Change this to the appropriate direcetory on your computer!!!


folderPath <- "Data/"
filename = "resampled Data"
resultsfilename = "results Data"
dir.create(filename)
dir.create(resultsfilename)
folderOut <- filename


studies <- dir(paste(filename, "/out", sep=""))

## only csv files
studies <- studies[grep(".csv", studies)]

files <- list.files(folderPath)
files <- files[grep(".csv", files)]
aucout <- matrix(NA, nrow=length(files), ncol=2)
for (j in 1:length(files)){
  
    source("./modeling/GWAS_functions.0,1,2genotype (attempt 4).R")
    genereteBatch(folderPath, folderOut, files[[j]])
    studies <- dir(paste(filename, "/out", sep=""))
    ## only csv files
    studies <- studies[grep(".csv", studies)]
    source("./modeling/_modelingFunctions_lda.R")
    batchAnalysis(paste(filename, "/out", sep=""), resultsfilename, studies[[j]], verbose = FALSE, outerFolds = 3L,
              outerRep = 2L)#, classifier = "classif.logreg") 
    studyname = substring(files[[j]], 1, 13)
    aucout[j,] <- c(studyname, aucanalysis)
    colnames(aucout) <- c("Study", "Simulated AUC")
}

View(aucout)


print("Done")

#--- end of the code
