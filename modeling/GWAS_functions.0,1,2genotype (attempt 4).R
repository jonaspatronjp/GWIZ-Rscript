# =====================================================#
# @author: Arnau Serra-Cayuela
# @contact: serracay@ualberta.ca
# @group: WishartLab (University of Alberta)
# @date: 2016
# @notes: code adapted partly from Joseph Landy (python)
# 
#
# From GWAS webpage information, the code generates
# a binary response depending on the frequence of 
# the allele and the odds ratio.
#
# ====================================================#


#--- set the working directory (select your own)
#setwd("/media/data/molecularU/GWAS_data_analysis/GWAS_scripts/dataIO/")
#setwd("/media/workspace/MolecularYouWKSpace/GWAS/GWAS_data_analysis/GWAS_scripts/dataIO/")

#setwd("/Users/wishart/Desktop/GWAS_project/")

# ================================================#
# =========== create disaggregated table =========

# Legend -----------------------------------#
# H --> Healthy (control, no exposure)
# D --> Disease (case, expoused)
# s --> Has the SNP polymorphism
# n --> No polymorphism
# ------------------------------------------#

### --- Function to estimate Hs, Hn, Ds and Dn from
#       the allele's frequency and the Odds Ratio.
#       This function is called from generateBinary
doTable <- function(freqPoly, Ncontrol, Ncase, OR , model){
  Hs <- ceiling(2 * Ncontrol * freqPoly) 
  # 'Healthy' (control / no exposure) non-polymorphism (Hn)
  Hn <- 2 * Ncontrol - Hs 
  # 'Disease' (case / exposure) with polymorphism (Ds)
  Ds <- ceiling( (OR * 2 * Ncase * Hs) / (Hn + OR * Hs) ) # developoing OR =(Ds*Hn)/(Dn*Hs) and Ncase = Ds + DN
  # 'Disesas' (case / exposure) non-polymorphism (Dn)
  Dn = 2 * Ncase - Ds
  # calc OR
  hatOR = ( (Ds*Hn)/(Dn*Hs) ) # estimated Odds Ratio
  return (list(Hs=Hs, Hn=Hn, Ds=Ds, Dn=Dn, hatOR=hatOR, model=model))
}

# example
#ORtable <- doTable(da[1, "risk_allele_freq"], da[1,"control.size"], da[1,"case.size"], da[1,"OR"])
#ORtable <- doTable(da[1, "risk.allele.freq"], da[1,"control.size"], da[1,"case.size"], da[1,"OR"])


### --- Function to check if the doTable created proper values
# --- A quick check for correctness:
# if cases/controls don't sum to case/control total, make a fuss warning.
# if OR do not coincide with the estimated OR, show a warning
# else: calculus has been done correctly
checkTable <- function(da, ORtable, SNP){
  #   if ( (ORtable$Hs + ORtable$Hn != da[1,"control.size"]) | (ORtable$Ds + ORtable$Dn != da[1,"case.size"]) ) {
  #     print ('Error: Hs + Hn != control size or Ds + Dn != case size ')
  #     proportion = FALSE
  #   } else{
  #     print ('Correct proportion')
  #     proportion = TRUE
  #   }
  if ( round(ORtable$hatOR, 2) == round(da[SNP,"OR"], 2) ){
    print ('Correct OR')
    ORcheck = TRUE
  } else{
    print ('Estimated odds ratio do not match with known OR')
    ORcheck = FALSE
  }
  #   if (proportion == TRUE & ORcheck == TRUE){
  #     invisible (TRUE)
  #   }else{
  #     invisible (FALSE)
  #   }
}

# example
#checkTable(da, ORtable)



# ================================================#
# =========== Generate Binary data =========



### --- Function to generate a binary output depending
#       on the values of Hs, Hn, Ds and Dn. Basically
#       it generates a vector of 1 and 0 depending if
#       there is a SNP or not, for the cases and the 
#       controls. Every vector is shuffled in order to
#       create the "observations"

doBinary <- function(ORtable){
  # create a vector of 1 (Ds) and 0 (TotalCases - Ds or Dn)
 
  # randomly shuffle the values
  set.seed(1984)# for reproductibility
  CasesbinaryVector <- c(rep(1, ORtable$Ds), rep(0, ORtable$Dn))
  CasesbinaryVector <- sample(CasesbinaryVector, length(CasesbinaryVector), replace = FALSE)
  #set.seed(1985)# for reproductibility
  CasesbinaryVector1 <- CasesbinaryVector[1:(length(CasesbinaryVector)/2)]
  CasesbinaryVector2 <- CasesbinaryVector[((length(CasesbinaryVector)/2)+1):length(CasesbinaryVector)]
  # for (i in 1: length(CasesbinaryVector)){
  #   if (sample(0:1, 1)==0){
  #     CasesbinaryVector[i] <- sample(0:1, 1)
  #   }
  # }
  # same for the controls
  # create a vector of 1 (Hs) and 0 (TotalControl - Hs or Hn)
  
  # for (i in 1: length(ControlbinaryVector)){
  #   if (sample(0:1, 1)==0){
  #     ControlbinaryVector[i] <- sample(0:1, 1)
  #   }
  # }
  # randomly shuffle the values
  set.seed(1984)# for reproductibility
  ControlbinaryVector <- c(rep(1, ORtable$Hs), rep(0, ORtable$Hn))
  ControlbinaryVector <- sample(ControlbinaryVector, length(ControlbinaryVector), replace = FALSE)
  #set.seed(1985)# for reproductibility
  ControlbinaryVector1 <- ControlbinaryVector[1:(length(ControlbinaryVector)/2)]
  ControlbinaryVector2 <- ControlbinaryVector[((length(ControlbinaryVector)/2)+1):length(ControlbinaryVector)]
  # merge cases and controls to one vector
  binaryVector1 <- c(CasesbinaryVector1, ControlbinaryVector1)
  binaryVector2 <- c(CasesbinaryVector2, ControlbinaryVector2)
  binaryVector <- c(rep(1, length(binaryVector1)))
  binaryVectorr <- c(rep(1, length(binaryVector1)))
  if (ORtable$model == "recessive"){
    print("recessive model")
  for (i in 1:length(binaryVector1)){
    if(as.numeric(binaryVector1[i])+as.numeric(binaryVector2[i])==2){
      binaryVector[i] <- 1
    }
    else{
    binaryVector[i] <- 0
    }}
  }
  else if (ORtable$model == "dominant"){
    print("dominant model")
  for (i in 1:length(binaryVector1)){
    if(as.numeric(binaryVector1[i])+as.numeric(binaryVector2[i])==0){
      binaryVector[i] <- 0
    }
    else{
      binaryVector[i] <- 1 
    }}
  }
  for (i in 1:length(binaryVector1)){
      binaryVectorr[i] <- as.numeric(binaryVector1[i])+as.numeric(binaryVector2[i])
    }
  return (binaryVector)
}
# x=0
# for (i in 1:length(ControlbinaryVector2)){
#   if(CasesbinaryVector1[i]==1){
#     x=x+1
#   }
# }
# print(x)
#example
#doBinary(ORtable)


### --- Function that outputs two dataframe:
#       1) Table with the Hs, Hn, Ds, Dn and predicted OR
#       2) Table with the 're-sampled' values

generateBinary <- function(dades, check = FALSE, filename){
  ## --- extract data (IMPORTANT: csv file MUST be consistent with the colum names)
  ## Number of control cases
  Ncontrol <- dades[1, 'control_size'] # usually it is located in column 8, 10(why double?)
  ## Number of case cases
  Ncase <- dades[1, 'case_size'] # usually it is located in column 7, 9(why double?)

  disaggregate <- list()
  binaries <- list()
  for (SNP in 1:dim(dades)[1] ){ # for each SNP
    
    ## since some files has the label risk_allele_freq and others risk.allele.freq
    if(sum(names(dades)=="risk_allele_freq")==0){ 
      names(dades)[which(names(dades)=="risk.allele.freq")] <- "risk_allele_freq"
    }
    print(dades)
    OR <- as.numeric(dades[SNP, "OR"])
    RAF <- as.numeric(dades[SNP, "risk_allele_freq"])
    
    if(OR<1){
      OR <- 1/OR
      RAF <- 1-RAF
    }
    
    #Ncontrol = dades[SNP, "control_size"]
    #Ncase = dades[SNP, "case_size"]
    model = dades[SNP, "model"]
    
    ORtable <- doTable (RAF, Ncontrol, Ncase, OR, model ) # dades[SNP, "risk.allele.freq"]
    
    disaggregate[[SNP]] <- ORtable
    if (check == TRUE){
      checkTable(dades, ORtable, SNP)
    }
    print(ORtable)
    binaryVector <- doBinary(ORtable)
    #binaryVector <- sample(binaryVector, length(binaryVector), replace = FALSE)
    # for (i in 1: length(binaryVector)){
    #   if (binaryVector[i] == 1 & sample(0:4, 1)==0){
    #     binaryVector[i] <- sample(0:1, 1)
    #   }
    # }
    # for (i in 1: length(binaryVector)){
    #   if (binaryVector[i] == 0 & sample(0:2, 1)==0){
    #     binaryVector[i] <- sample(0:1, 1)
    #   }
    # }
    binaries[[SNP]] <- binaryVector
  }
  
  ## format disaggregate to output a data frame with the Hs, Hn... values
  disaggregate <- as.data.frame(t(sapply(disaggregate, function (x) unlist(x)))) # data.frame 1)
  row.names(disaggregate) <- dades$accession
  
  # data frame from the collected binaries list object
  binariesDF <- as.data.frame(binaries)
  
  #--- variable name row (accession number + phenotype + GWAScode)  # GWAS is obtained from the filename
  vars <- paste( dades[, "accession"]) #, dades[1, "phenotype"], gsub("\\..*", "", filename), sep = '_')
  names(binariesDF) <- vars
  #--- create individuals label (observation)
  subject <- paste( "s", seq(1, Ncontrol + Ncase), sep = "" )
  #--- create class type label (case, control...)
  group <- as.factor( c( rep("case", Ncase), rep("control", Ncontrol) ) )
  
  df <- data.frame(subject = subject, class = group, binariesDF) # data frame with the binary values. data.frame 2)
  return (list(disaggregate, df))
}

# example
#generateBinary(da, filename='whateverstudyis')


#---  Function to perform the analysis in batch. Only a folder path is required.
#     That folfer may content the csv files that are properly formated (have the appropriate columns)
genereteBatch <- function(folderPath, folderOut, files){
  ## create output directory if it does not already exist
  if (!dir.exists(file.path(folderOut, "out"))){
    dir.create(file.path(folderOut, "out"), showWarnings = FALSE)
    print ("Folder 'out' created")
  }
  if (!dir.exists(file.path(folderOut, "outDisaggregate"))){
    dir.create(file.path(folderOut, "outDisaggregate"), showWarnings = FALSE)
    print ("Folder 'outDisaggregate' created")
  }
  
  #files <- list.files(folderPath)
  #files <- files[-which(files=="out")]
  #files <- files[-which(files=="outDisaggregate")]
  #files <- files[grep(".csv", files)]

  #files = c("Combining Information from Common Type 2 Diabetes Risk Polymorphisms Improves Disease Prediction.csv")
  cat ("Resampling ", length(files), " studies:\n")
  print (files)
  print ('=================================================================================')
  for (study in files){
    da <- read.csv(file.path(folderPath, study))
    da[is.na(da)]<-""
    print (study)
    out <- generateBinary(da, filename = study, check = FALSE)
    
    # set the files names
    #diseaseName <- da$phenotype[1]
    #GWAScode <- gsub("\\..*", "", study)
    #fileNameBase <- paste(diseaseName, GWAScode, sep = "_")
    
    fileNameBase <- gsub(".csv", "", study)
    
    #fileName <- paste(gsub("\\..*", "", study), "sampling", sep = "_")
    #diseaseName <- paste(da$phenotype,fileName, sep= "")
    fileName1 <- paste(fileNameBase, "TableDisaggregate.csv", sep = "_")
    fileName2 <- paste(fileNameBase, "Resampled.csv", sep = "_")
    # export the data with the table of disaggregated values (Hs, Hn...)
    write.csv(data.frame(da, out[[1]]), file.path(folderOut, "outDisaggregate", fileName1), row.names = FALSE)
    # export the resamped sample (binary values)
    write.csv(out[[2]], file.path(folderOut, "out", fileName2), row.names = FALSE)
  }
  print("done!!! :)")
}

# Example
#genereteBatch(folderpath)
# It may output some warning messages, which are associates with that some of the names have spaces or special characters (ex:barret's), but performs correctly the analysis







#### --- end of the code

######################################################################################################
######################################################################################################
### ================================   DEBBUG 
######################################################################################################
######################################################################################################

# 
# # --- extract data (IMPORTANT: csv file MUST be consistent with the colum names)
# # Number of control cases
# Ncontrol <- da[1, 'control.size'] # usually it is located in column 8, 10(why double?)
# # Number of case cases
# Ncase <- da[1, 'case.size'] # usually it is located in column 7, 9(why double?)
# 
# 
# # Legend -----------------------------------#
# # H --> Healthy (control, no exposure)
# # D --> Disease (case, expoused)
# # s --> Has the SNP polymorphism
# # n --> No polymorphism
# # ------------------------------------------#
# 
# # 'Healthy' (control / no exposure) with polymorphism:(Hs)
# freqPoly = da[1, "risk_allele_freq"] # frequency of the polymorphism
# Hs <- ceiling(Ncontrol * freqPoly) ; Hs
# 
# # 'Healthy' (control / no exposure) non-polymorphism (Hn)
# Hn <- Ncontrol - Hs ; Hs
# 
# # 'Disease' (case / exposure) with polymorphism (Ds)
# OR<-da[1, "OR"] # Odds ratio
# Ds <- ceiling( (OR * Ncase * Hs) / (Hn + OR*Hs) ) # developoing OR =(Ds*Hn)/(Dn*Hs) and Ncase = Ds + DN
# 
# # 'Disease' (case / exposure) non-polymorphism (Dn)
# Dn = Ncase - Ds
# 
# # 're'-calculate OR, to see if match...
# ORcheck = (Ds*Hn)/(Dn*Hs)
# round(OR,2) == round(ORcheck,2)
# Hn/Hs
# 
# Hs; Hn; Ds; Dn
# 
# 
# 
# # create a vector of 1 (Ds) and 0 (TotalCases - Ds or Dn)
# CasesbinaryVector <- c(rep(1, ORtable$Ds),rep(0, ORtable$Dn))
# # randomly shuffle the values
# set.seed(1984)# for reproductibility
# CasesbinaryVector <- sample(CasesbinaryVector, length(CasesbinaryVector), replace = FALSE)
# # same for the controls
# # create a vector of 1 (Hs) and 0 (TotalControl - Hs or Hn)
# ControlbinaryVector <- c(rep(1, ORtable$Hs),rep(0, ORtable$Hn))
# # randomly shuffle the values
# set.seed(1984)# for reproductibility
# ControlbinaryVector <- sample(ControlbinaryVector, length(ControlbinaryVector), replace = FALSE)
# # merge cases and controls to one vector
# binaryVector <- c(CasesbinaryVector, ControlbinaryVector)
# 
# 
# #--- create individuals label (observation)
# subject <- paste( "s", seq(1, Ncontrol + Ncase), sep = "" )
# 
# #--- create class type label (case, control...)
# group <- as.factor( c( rep("case", Ncase), rep("control", Ncontrol) ) )
# 
# #--- variable name row (accession number + phenotype + GWAScode)  # GWAS is obtained from the filename
# vars <- paste( da[, "accession"], da[1, "phenotype"], gsub("\\..*", "", filename), sep = '_')
# 
# 
# 
# 
# 
# 
# # create a data frame
# result <- data.frame( subject = subject, class = group, binaryVector)
# 
# 
# 
# #### ------ all in one
# 
# check = FALSE
# 
# 
# disaggregate <- list()
# binaries <- list()
# for (SNP in 1:dim(da)[1] ){
#   print (SNP)
#   ORtable <- doTable (da[SNP, "risk_allele_freq"], Ncontrol, Ncase, da[SNP, "OR"] )
#   disaggregate[[SNP]] <- ORtable
#   if (check == TRUE){
#     checkTable(ORtable, OR = da[SNP, "OR"])
#   }
#   binaryVector <- doBinary(ORtable)
#   binaries[[SNP]] <- binaryVector
# }
# binariesDF <- as.data.frame(binaries)
# #--- variable name row (accession number + phenotype + GWAScode)  # GWAS is obtained from the filename
# vars <- paste( da[, "accession"], da[1, "phenotype"], gsub("\\..*", "", filename), sep = '_')
# names(binariesDF) <- vars
# #--- create individuals label (observation)
# subject <- paste( "s", seq(1, Ncontrol + Ncase), sep = "" )
# #--- create class type label (case, control...)
# group <- as.factor( c( rep("case", Ncase), rep("control", Ncontrol) ) )
# 
# df <- data.frame(subject = subject, class = group, binariesDF)
# 
# # -- as function




##############################################################
# =============  Debugging

# 
# folderpath <- "/media/data/molecularU/GWAS_data_analysis/GWAS_scripts/dataIO/Studies"
# files <- list.files(Filespath)
# files <- files[-which(files=="out")]
# 
# for (file in files){
#   print (file)
#   da <- read.csv(file.path("Studies",file))
#   out <- generateBinary(da, filename = file)
#   fileName <- paste(gsub("\\..*", "", file), "sampling", sep = "_")
#   fileName <- paste(fileName,"csv", sep = ".")
#   #write.csv(out, file.path("./..",Filespath, "out", fileName))
#   write.csv(out, file.path(folderpath, "out", fileName), row.names = FALSE)
# }
# print("done")
# 
# 
# 
# 
# if (!dir.exists(file.path(folderpath, "out"))){
#   print ("si")
# }
# 
# 
# 
# study <- "HGVST163.csv"
# da <- read.csv(file.path(folderpath, study))
# out <- generateBinary(da, filename = study, check = FALSE)
# fileName <- paste(gsub("\\..*", "", study), "sampling", sep = "_")
# fileName <- paste(fileName, "csv", sep = ".")
# #write.csv(out, file.path("./..",Filespath, "out", fileName))
# write.csv(out[[2]], file.path(folderpath, "out", fileName), row.names = FALSE)
# 

