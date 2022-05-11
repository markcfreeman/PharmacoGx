rm(list = ls())
cat("\f")

options(error = recover)
setwd("~/Skainki")

library(nleqslv)
library(Rcpp)

source("drugCellLineCombinationAnalysis.R")
source("katzirModel.R")
sourceCpp("C:\\Users\\Marky\\Documents\\Skainki\\Katzir\\katzirMeshOptimization.cpp")

tissueType <- "lung"
dataset <- "MERCK"
modelType <- "unidirectional"
maxResidual <- 0.2
maxObservedViabilityFloor <- 0.9
candidateTripletViabilityCeiling <- 0.1

load(paste0(tissueType, "_", dataset, "_singleDrugData.RData"))
load(paste0(tissueType, "_", dataset, "_twoDrugData.RData"))
load(paste0(tissueType, "_", dataset, "_usableDrugCellLineCombinations.RData"))

# comparisonTableColNames <- c("Cell Line",
#                              "Drug A",
#                              "Drug B",
#                              "logD0i",
#                              "logD0j",
#                              "logNi",
#                              "logNj",
#                              "logAij",
#                              "logAji",
#                              "residual")
#comparisonTable <- matrix(NA, nrow = 3 * combosToSample, ncol = length(comparisonTableColNames))
#colnames(comparisonTable) <- comparisonTableColNames
pairInd <- 1
tripletInd <- 1
candidateTriplets <- as.list(rep(NA, length(usableDrugCellLineCombinations)))
candidatePairs <- candidateTriplets

for (i in seq_along(usableDrugCellLineCombinations))
{
  combination <- usableDrugCellLineCombinations[[i]]
  combinationAnalysis <- analyzeInteraction(drugs = unlist(combination@drugs),
                                            cellLine = combination@cellLine,
                                            tissueType = tissueType,
                                            dataset = dataset,
                                            modelType = modelType)
  
  print(c(unlist(combination@drugs), combination@cellLine))
  isCandidateTriplet <- TRUE
  
  for (A in 1:(length(unlist(combination@drugs)) - 1))
  {
    for (B in (A + 1):length(unlist(combination@drugs)))
    {
      drugA <- unlist(combination@drugs)[A]
      drugB <- unlist(combination@drugs)[B]
      drugPair <- paste(drugA, "&", drugB)
      
      # if (combinationAnalysis$katzirResiduals[[drugPair]] > maxResidual)
      # {
      #   isCandidateTriplet <- FALSE
      #   print(paste(drugPair,
      #         ": excessive residual (",
      #         combinationAnalysis$katzirResiduals[[drugPair]],
      #         ")"))
      # }
      
      if ("logAij" %in% names(combinationAnalysis$katzirModelFits[[drugPair]]) &&
               combinationAnalysis$katzirModelFits[[drugPair]][["logAij"]] > 0)
      {
        print(paste(drugPair,
                    "antagonistic interaction (",
                    combinationAnalysis$katzirModelFits[[drugPair]][["logAij"]],
                    ")"))
        isCandidateTriplet <- FALSE
        next
      }
      
      if ("logAji" %in% names(combinationAnalysis$katzirModelFits[[drugPair]]) &&
               combinationAnalysis$katzirModelFits[[drugPair]][["logAji"]] > 0)
      {
        print(paste(drugPair,
              "antagonistic interaction (",
              combinationAnalysis$katzirModelFits[[drugPair]][["logAji"]],
              ")"))
        isCandidateTriplet <- FALSE
        next
      }
      
      # for (drug in c(drugA, drugB))
      # {
      #   singleDrugDoseResponseData <- extractDrugCellLineCombinationData(drugs = list(drug),
      #                                                                    cellLine = combination@cellLine,
      #                                                                    allData = singleDrugData,
      #                                                                    dataType = "singleDrug")
      #   
      #   if (min(singleDrugDoseResponseData[, "Viability"]) > maxObservedViabilityFloor)
      #   {
      #     isCandidateTriplet <- FALSE
      #     print(paste(drug,
      #           ": insufficient observed response (",
      #           min(singleDrugDoseResponseData[, "Viability"]),
      #           ")"))
      #   }
      # }
      
      data <- extractDrugCellLineCombinationData(drugs = as.list(c(drugA, drugB)),
                                                 cellLine = combination@cellLine,
                                                 allData = twoDrugData,
                                                 dataType = "twoDrug")
      CIs <- bootstrapParamCIs(data = data,
                               modelType = modelType,
                               fittedParams = combinationAnalysis$katzirModelFits[[drugPair]],
                               hillModelFits = list(combinationAnalysis$hillModelFits[[drugA]],
                                                    combinationAnalysis$hillModelFits[[drugB]]),
                               hillModelFitConfidenceIntervals = list(combinationAnalysis$hillModelFitConfidenceIntervals[[drugA]],
                                                                      combinationAnalysis$hillModelFitConfidenceIntervals[[drugB]]))
      
      if ("logAij" %in% names(combinationAnalysis$katzirModelFits[[drugPair]]) &&
          CIs["logAij", "97.5%"] > 0)
      {
        print(paste(drugA, "&", drugB, "uncertain interaction character"))
        isCandidateTriplet <- FALSE
        next
      }
      
      if ("logAji" %in% names(combinationAnalysis$katzirModelFits[[drugPair]]) &&
          CIs["logAji", "97.5%"] > 0)
      {
        print(paste(drugA, "&", drugB, "uncertain interaction character"))
        isCandidateTriplet <- FALSE
        next
      }
      
      print(paste(drugA, "&", drugB, "synergy"))
      candidatePairs[[pairInd]] <- list(drugA, drugB, combination@cellLine)
      pairInd <- pairInd + 1
    }
  }
  
  if (isCandidateTriplet)
  {
    concentrations <- list()
    
    for (drug in unlist(combination@drugs))
    {
      singleDrugDoseResponseData <- extractDrugCellLineCombinationData(drugs = list(drug),
                                                                       cellLine = combination@cellLine,
                                                                       allData = singleDrugData,
                                                                       dataType = "singleDrug")
      concentrations[[drug]] <- unique(singleDrugDoseResponseData[, "Concentration"])
    }
    
    x <- expand.grid(DrugAConcentration = concentrations[[1]],
                     DrugBConcentration = concentrations[[2]],
                     DrugCConcentration = concentrations[[3]],
                     stringsAsFactors = FALSE)
    logD0s <- c(combinationAnalysis$hillModelFits[[1]]$m$getAllPars()[["logD0"]],
                combinationAnalysis$hillModelFits[[2]]$m$getAllPars()[["logD0"]],
                combinationAnalysis$hillModelFits[[3]]$m$getAllPars()[["logD0"]])
    logNs <- c(combinationAnalysis$hillModelFits[[1]]$m$getAllPars()[["logN"]],
               combinationAnalysis$hillModelFits[[2]]$m$getAllPars()[["logN"]],
               combinationAnalysis$hillModelFits[[3]]$m$getAllPars()[["logN"]])
    logAs <- matrix(0, nrow = 3, ncol = 3)
    
    if ("logAij" %in% names(combinationAnalysis$katzirModelFits[[1]]))
    {
      logAs[1, 2] <- combinationAnalysis$katzirModelFits[[1]][["logAij"]]
    }
    
    if ("logAij" %in% names(combinationAnalysis$katzirModelFits[[2]]))
    {
      logAs[1, 3] <- combinationAnalysis$katzirModelFits[[2]][["logAij"]]
    }
    
    if ("logAji" %in% names(combinationAnalysis$katzirModelFits[[1]]))
    {
      logAs[2, 1] <- combinationAnalysis$katzirModelFits[[1]][["logAji"]]
    }
    
    if ("logAij" %in% names(combinationAnalysis$katzirModelFits[[3]]))
    {
      logAs[2, 3] <- combinationAnalysis$katzirModelFits[[3]][["logAij"]]
    }
    
    if ("logAji" %in% names(combinationAnalysis$katzirModelFits[[2]]))
    {
      logAs[3, 1] <- combinationAnalysis$katzirModelFits[[2]][["logAji"]]
    }
    
    if ("logAji" %in% names(combinationAnalysis$katzirModelFits[[3]]))
    {
      logAs[3, 2] <- combinationAnalysis$katzirModelFits[[3]][["logAji"]]
    }
    
    viabilities <- lapply(seq_len(nrow(x)),
                          function(i){ return(KatzirUnidirectional(x = x[i, ],
                                                                   logD0s = logD0s,
                                                                   logNs = logNs,
                                                                   logAs = logAs)) })
    
    if (any(unlist(viabilities) < candidateTripletViabilityCeiling))
    {
      print("****************************************")
      print("Potential candidate triplet identified!")
      #print(combinationAnalysis$katzirModelFits)
      print("****************************************")
      
      candidateTriplets[[tripletInd]] <- combinationAnalysis
      tripletInd <- tripletInd + 1
    }
  }
  
  cat("\n")
}