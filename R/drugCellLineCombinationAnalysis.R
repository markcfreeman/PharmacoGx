extractDrugCellLineCombinationData <- function(drugs, cellLine, allData, dataType)
{
  if (dataType == "singleDrug")
  {
    drugDataRows <- which(allData[, "Drug"] %in% drugs & allData[, "Cell Line"] == cellLine)
    relevantColNames <- c("Concentration", "Viability")
  }
  else if (dataType == "twoDrug")
  {
    drugDataRows <- which(allData[, "Drug A"] %in% drugs &
                            allData[, "Drug B"] %in% drugs &
                            allData[, "Cell Line"] == cellLine)
    relevantColNames <- c("Drug A Concentration", "Drug B Concentration", "Viability")
  }
  else
  {
    stop("Invalid dataType.")
  }
  
  drugDataCols <- which(colnames(allData) %in% relevantColNames)
  doseResponseData <- allData[drugDataRows, drugDataCols]
  colnames(doseResponseData) <- relevantColNames
  
  return(doseResponseData)
}

multiDrugCellLineCombinationAnalysis <- setRefClass("multiDrugCellLineCombinationAnalysis",
                                                    fields = list(drugs = "list",
                                                                  cellLine = "character",
                                                                  tissueType = "character",
                                                                  dataset = "character",
                                                                  hillModelFits = "list",
                                                                  katzirModelFits = "list"),
                                                    methods = list(
                                                      drugPairs = function()
                                                      {
                                                        drugPairs <- list()
                                                        
                                                        for (i in 1:(length(drugs) - 1))
                                                        {
                                                          for (j in (i + 1):length(drugs))
                                                          {
                                                            drugPairs <- c(drugPairs, paste(drugs[[i]], "&", drugs[[j]]))
                                                          }
                                                        }
                                                        
                                                        return(drugPairs)
                                                      },
                                                      plotResults = function()
                                                      {
                                                        for (drug in drugs)
                                                        {
                                                          xrange <- range(singleDrugDoseResponseData[[drug]][, "Concentration"])
                                                          interpolatedLogConcentrations <- seq(from = log(min(xrange)), to = log(max(xrange)), length.out = 1001)
                                                          plot(NULL,
                                                               log = "x",
                                                               xlim = xrange,
                                                               ylim = c(0, 1),
                                                               main = paste(drug, ":", cellLine),
                                                               xlab = "Concentration (uM)",
                                                               ylab = "Viability")
                                                          points(singleDrugDoseResponseData[[drug]], pch = 21, col = "blue")
                                                          points(exp(interpolatedLogConcentrations), Hill(exp(interpolatedLogConcentrations),
                                                                                                          hillModelFits[[drug]]$m$getPars()[["logD0"]],
                                                                                                          hillModelFits[[drug]]$m$getPars()[["logN"]]),
                                                                 pch = 20,
                                                                 cex = 0.1,
                                                                 col = "blue")
                                                          abline(h = 0.1,
                                                                 pch = 20,
                                                                 col = "red")
                                                        }
                                                        
                                                        for (drugPair in drugPairs())
                                                        {
                                                          xs <- sort(unique(twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"]))
                                                          ys <- sort(unique(twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"]))
                                                          logXRange <- log(range(xs))
                                                          logYRange <- log(range(ys))
                                                          interpolatedLogXs <- seq(from = logXRange[1] - 0.5, to = logXRange[2] + 0.5, length.out = 101)
                                                          interpolatedLogYs <- seq(from = logYRange[1] - 0.5, to = logYRange[2] + 0.5, length.out = 101)
                                                          observedViabilities <- matrix(twoDrugDoseResponseData[[drugPair]][, "Viability"], nrow = 4, byrow = TRUE)
                                                          predictedViabilities <- outer(interpolatedLogXs,
                                                                                        interpolatedLogYs,
                                                                                        function(x, y){return(Katzir2(exp(x),
                                                                                                                      exp(y),
                                                                                                                      katzirModelFits[[drugPair]]$m$getPars()[["logD0i"]],
                                                                                                                      katzirModelFits[[drugPair]]$m$getPars()[["logD0j"]],
                                                                                                                      katzirModelFits[[drugPair]]$m$getPars()[["logNi"]],
                                                                                                                      katzirModelFits[[drugPair]]$m$getPars()[["logNj"]],
                                                                                                                      katzirModelFits[[drugPair]]$m$getPars()[["logAij"]],
                                                                                                                      katzirModelFits[[drugPair]]$m$getPars()[["logAji"]]))})
                                                          colourScale <- rev(rainbow(100))
                                                          image(x = interpolatedLogXs,
                                                                y = interpolatedLogYs,
                                                                z = predictedViabilities,
                                                                zlim = c(0, 1),
                                                                main = paste(drugPair, ":", cellLine),
                                                                xlab = "Drug A Log Concentration (uM)",
                                                                ylab = "Drug B Log Concentration (uM)",
                                                                col = colourScale)
                                                          points(x = log(twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"]),
                                                                 y = log(twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"]),
                                                                 pch = 21,
                                                                 col = "black",
                                                                 bg = colourScale[ceiling(twoDrugDoseResponseData[[drugPair]][, "Viability"] *
                                                                                            length(colourScale))])
                                                        }
                                                      }))

fitHillModelToData <- function(data, maxLogN = 2, minLogN = -maxLogN)
{
  minLogD0 <- log(min(data[, "Concentration"])) - 1
  maxLogD0 <- log(max(data[, "Concentration"])) + 1
  hillModelFit <- nls(formula = y ~ Hill(x, logD0, logN),
                      data = list(y = data[, "Viability"],
                                  x = data[, "Concentration"]),
                      start = list(logD0 = log(data[which.min(abs(data[, "Viability"] - 0.5)), "Concentration"]),
                                   logN = 0),
                      lower = list(logD0 = minLogD0, logN = minLogN),
                      upper = list(logD0 = maxLogD0, logN = maxLogN),
                      algorithm = "port",
                      control = nls.control(warnOnly = TRUE))
  hillModelFitConfidenceIntervals <- matrix(c(-Inf, minLogN, Inf, maxLogN),
                                            nrow = 2,
                                            dimnames = list(c("logD0", "logN"), c("2.5%", "97.5%")))
  try(expr = {hillModelFitConfidenceIntervals <- confint(hillModelFit)}, TRUE)
  
  if (is.na(hillModelFitConfidenceIntervals["logD0", "2.5%"]))
  {
    hillModelFitConfidenceIntervals["logD0", "2.5%"] <- minLogD0
  }
  
  if (is.na(hillModelFitConfidenceIntervals["logD0", "97.5%"]))
  {
    hillModelFitConfidenceIntervals["logD0", "97.5%"] <- maxLogD0
  }
  
  if (is.na(hillModelFitConfidenceIntervals["logN", "2.5%"]))
  {
    hillModelFitConfidenceIntervals["logN", "2.5%"] <- minLogN
  }
  
  if (is.na(hillModelFitConfidenceIntervals["logN", "97.5%"]))
  {
    hillModelFitConfidenceIntervals["logN", "97.5%"] <- maxLogN
  }
  
  return(list(data = data, fit = hillModelFit, CIs = hillModelFitConfidenceIntervals))
}

fitKatzirModelToData <- function(data,
                                 hillModelFits,
                                 hillModelFitConfidenceIntervals,
                                 modelType = "unidirectional",
                                 maxLogN = 2,
                                 minLogN = -maxLogN,
                                 maxLogA = 5,
                                 minLogA = -maxLogA)
{
  if (modelType == "lowData")
  {
    minLogD0i <- log(min(data[, "Drug A Concentration"])) - 1
    minLogD0j <- log(min(data[, "Drug B Concentration"])) - 1
    maxLogD0i <- log(max(data[, "Drug A Concentration"])) + 1
    maxLogD0j <- log(max(data[, "Drug B Concentration"])) + 1
    
    start1 <- list(logD0i = log(data[which.min(abs(data[, "Viability"] - 0.25)), "Drug A Concentration"]),
                   logD0j = log(data[which.min(abs(data[, "Viability"] - 0.25)), "Drug B Concentration"]),
                   logAji = 0)
    lower1 <- list(logD0i = minLogD0i,
                   logD0j = minLogD0j,
                   logAji = minLogA)
    upper1 <- list(logD0i = maxLogD0i,
                   logD0j = maxLogD0j,
                   logAji = maxLogA)
    start1 <- refineKatzirStart(y = data[, "Viability"],
                                xi = data[, "Drug A Concentration"],
                                xj = data[, "Drug B Concentration"],
                                start = start1,
                                lower = lower1,
                                upper = upper1,
                                missingParams = c("logNi", "logNj", "logAij"))
    startResidual1 <- sum((y - Katzir2(xi = data[, "Drug A Concentration"],
                                       xj = data[, "Drug B Concentration"],
                                       logD0i = start1[["logD0i"]],
                                       logD0j = start1[["logD0j"]],
                                       logNi = 0,
                                       logNj = 0,
                                       logAij = 0,
                                       logAji = start1[["logAji"]])) ^ 2)
    fitResidual1 <- Inf
    try(expr = {fit1 <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, 0, 0, 0, logAji),
                            data = list(y = data[, "Viability"],
                                        xi = data[, "Drug A Concentration"],
                                        xj = data[, "Drug B Concentration"]),
                            start = start1,
                            lower = lower1,
                            upper = upper1,
                            algorithm = "port",
                            control = nls.control(warnOnly = TRUE));
    fitResidual1 <- fit1$m$deviance()},
        TRUE)
    
    start2 <- list(logD0i = log(data[which.min(abs(data[, "Viability"] - 0.25)), "Drug A Concentration"]),
                   logD0j = log(data[which.min(abs(data[, "Viability"] - 0.25)), "Drug B Concentration"]),
                   logAji = 0)
    lower2 <- list(logD0i = minLogD0i,
                   logD0j = minLogD0j,
                   logAij = minLogA)
    upper2 <- list(logD0i = maxLogD0i,
                   logD0j = maxLogD0j,
                   logAij = maxLogA)
    start2 <- refineKatzirStart(y = data[, "Viability"],
                                xi = data[, "Drug A Concentration"],
                                xj = data[, "Drug B Concentration"],
                                start = start2,
                                lower = lower2,
                                upper = upper2,
                                missingParams = c("logNi", "logNj", "logAji"))
    startResidual2 <- sum((data[, "Viability"] -
                             Katzir2(xi = data[, "Drug A Concentration"],
                                     xj = data[, "Drug B Concentration"],
                                     logD0i = start2[["logD0i"]],
                                     logD0j = start2[["logD0j"]],
                                     logNi = 0,
                                     logNj = 0,
                                     logAij = start2[["logAij"]],
                                     logAji = 0)) ^ 2)
    fitResidual2 <- Inf
    try(expr = {fit2 <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, 0, 0, logAij, 0),
                            data = list(y = data[, "Viability"],
                                        xi = data[, "Drug A Concentration"],
                                        xj = data[, "Drug B Concentration"]),
                            start = start2,
                            lower = lower2,
                            upper = upper2,
                            algorithm = "port",
                            control = nls.control(warnOnly = TRUE));
    fitResidual2 <- fit2$m$deviance()},
        TRUE)
    
    minResidual <- min(c(startResidual1, startResidual2, fitResidual1, fitResidual2))
    
    if (startResidual1 == minResidual)
    {
      katzirModelFit <- start1
    }
    else if (startResidual2 == minResidual)
    {
      katzirModelFit <- start2
    }
    else if (fitResidual1 == minResidual)
    {
      katzirModelFit <- fit1$m$getPars() 
    }
    else
    {
      katzirModelFit <- fit2$m$getPars() 
    }
  }
  else if (modelType == "unidirectional")
  {
    start1 <- list(logD0i = hillModelFits[[1]]$m$getPars()[["logD0"]],
                   logD0j = hillModelFits[[2]]$m$getPars()[["logD0"]],
                   logNi = hillModelFits[[1]]$m$getPars()[["logN"]],
                   logNj = hillModelFits[[2]]$m$getPars()[["logN"]],
                   logAij = 0)
    lower1 <- list(logD0i = hillModelFitConfidenceIntervals[[1]]["logD0", "2.5%"],
                  logD0j = hillModelFitConfidenceIntervals[[2]]["logD0", "2.5%"],
                  logNi = hillModelFitConfidenceIntervals[[1]]["logN", "2.5%"],
                  logNj = hillModelFitConfidenceIntervals[[2]]["logN", "2.5%"],
                  logAij = minLogA)
    upper1 <- list(logD0i = hillModelFitConfidenceIntervals[[1]]["logD0", "97.5%"],
                  logD0j = hillModelFitConfidenceIntervals[[2]]["logD0", "97.5%"],
                  logNi = hillModelFitConfidenceIntervals[[1]]["logN", "97.5%"],
                  logNj = hillModelFitConfidenceIntervals[[2]]["logN", "97.5%"],
                  logAij = maxLogA)
    start1 <- refineKatzirStart(y = data[, "Viability"],
                                xi = data[, "Drug A Concentration"],
                                xj = data[, "Drug B Concentration"],
                                start = start1,
                                lower = lower1,
                                upper = upper1,
                                missingParams = c("logAji"))
    startResidual1 <- sum((data[, "Viability"] -
                             Katzir2(xi = data[, "Drug A Concentration"],
                                     xj = data[, "Drug B Concentration"],
                                     logD0i = start1[["logD0i"]],
                                     logD0j = start1[["logD0j"]],
                                     logNi = start1[["logNi"]],
                                     logNj = start1[["logNj"]],
                                     logAij = start1[["logAij"]],
                                     logAji = 0)) ^ 2)
    
    fitResidual1 <- Inf
    try(expr = {fit1 <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, logNi, logNj, logAij, 0),
                            data = list(y = data[, "Viability"],
                                        xi = data[, "Drug A Concentration"],
                                        xj = data[, "Drug B Concentration"]),
                            start = start1,
                            lower = lower1,
                            upper = upper1,
                            algorithm = "port");
    fitResidual1 <- fit1$m$deviance()},
        TRUE)
    
    start2 <- list(logD0i = hillModelFits[[1]]$m$getPars()[["logD0"]],
                   logD0j = hillModelFits[[2]]$m$getPars()[["logD0"]],
                   logNi = hillModelFits[[1]]$m$getPars()[["logN"]],
                   logNj = hillModelFits[[2]]$m$getPars()[["logN"]],
                   logAji = 0)
    lower2 <- list(logD0i = hillModelFitConfidenceIntervals[[1]]["logD0", "2.5%"],
                  logD0j = hillModelFitConfidenceIntervals[[2]]["logD0", "2.5%"],
                  logNi = hillModelFitConfidenceIntervals[[1]]["logN", "2.5%"],
                  logNj = hillModelFitConfidenceIntervals[[2]]["logN", "2.5%"],
                  logAji = minLogA)
    upper2 <- list(logD0i = hillModelFitConfidenceIntervals[[1]]["logD0", "97.5%"],
                  logD0j = hillModelFitConfidenceIntervals[[2]]["logD0", "97.5%"],
                  logNi = hillModelFitConfidenceIntervals[[1]]["logN", "97.5%"],
                  logNj = hillModelFitConfidenceIntervals[[2]]["logN", "97.5%"],
                  logAji = maxLogA)
    start2 <- refineKatzirStart(y = data[, "Viability"],
                                xi = data[, "Drug A Concentration"],
                                xj = data[, "Drug B Concentration"],
                                start = start2,
                                lower = lower2,
                                upper = upper2,
                                missingParams = c("logAij"))
    startResidual2 <- sum((data[, "Viability"] -
                             Katzir2(xi = data[, "Drug A Concentration"],
                                     xj = data[, "Drug B Concentration"],
                                     logD0i = start2[["logD0i"]],
                                     logD0j = start2[["logD0j"]],
                                     logNi = start2[["logNi"]],
                                     logNj = start2[["logNj"]],
                                     logAij = 0,
                                     logAji = start2[["logAji"]])) ^ 2)
    fitResidual2 <- Inf
    try(expr = {fit2 <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, logNi, logNj, 0, logAji),
                            data = list(y = data[, "Viability"],
                                        xi = data[, "Drug A Concentration"],
                                        xj = data[, "Drug B Concentration"]),
                            start = start2,
                            lower = lower2,
                            upper = upper2,
                            algorithm = "port");
    fitResidual2 <- fit2$m$deviance()},
        TRUE)
    
    minResidual <- min(c(fitResidual1, fitResidual2, startResidual1, startResidual2))
    
    if (startResidual1 == minResidual)
    {
      katzirModelFit <- start1
    }
    else if (startResidual2 == minResidual)
    {
      katzirModelFit <- start2
    }
    else if (fitResidual1 == minResidual)
    {
      katzirModelFit <- fit1$m$getPars()
    }
    else
    {
      katzirModelFit <- fit2$m$getPars()
    }
  }
  else if (modelType == "bidirectional")
  {
    start <- list(logD0i = hillModelFits[[1]]$m$getPars()[["logD0"]],
                  logD0j = hillModelFits[[2]]$m$getPars()[["logD0"]],
                  logNi = hillModelFits[[1]]$m$getPars()[["logN"]],
                  logNj = hillModelFits[[2]]$m$getPars()[["logN"]],
                  logAij = 0,
                  logAji = 0)
    lower <- list(logD0i = max(hillModelFitConfidenceIntervals[[1]]["logD0", "2.5%"],
                               minLogD0i,
                               na.rm = TRUE),
                  logD0j = max(hillModelFitConfidenceIntervals[[2]]["logD0", "2.5%"],
                               minLogD0j,
                               na.rm = TRUE),
                  logNi = max(hillModelFitConfidenceIntervals[[1]]["logN", "2.5%"],
                              minLogN,
                              na.rm = TRUE),
                  logNj = max(hillModelFitConfidenceIntervals[[2]]["logN", "2.5%"],
                              minLogN,
                              na.rm = TRUE),
                  logAij = minLogA,
                  logAji = minLogA)
    upper <- list(logD0i = min(hillModelFitConfidenceIntervals[[1]]["logD0", "97.5%"],
                               maxLogD0i,
                               na.rm = TRUE),
                  logD0j = min(hillModelFitConfidenceIntervals[[2]]["logD0", "97.5%"],
                               maxLogD0j,
                               na.rm = TRUE),
                  logNi = min(hillModelFitConfidenceIntervals[[1]]["logN", "97.5%"],
                              maxLogN,
                              na.rm = TRUE),
                  logNj = min(hillModelFitConfidenceIntervals[[2]]["logN", "97.5%"],
                              maxLogN,
                              na.rm = TRUE),
                  logAij = maxLogA,
                  logAji = maxLogA)
    start <- refineKatzirStart(y = data[, "Viability"],
                               xi = data[, "Drug A Concentration"],
                               xj = data[, "Drug B Concentration"],
                               start = start,
                               lower = lower,
                               upper = upper)
    startResidual <- sum((data[, "Viability"] -
                            Katzir2(data[, "Drug A Concentration"],
                                    data[, "Drug B Concentration"],
                                    start[["logD0i"]],
                                    start[["logD0j"]],
                                    start[["logNi"]],
                                    start[["logNj"]],
                                    start[["logAij"]],
                                    start[["logAji"]])) ^ 2)
    fitResidual <- Inf
    try(expr = {fit <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, logNi, logNj, logAij, logAji),
                           data = list(y = data[, "Viability"],
                                       xi = data[, "Drug A Concentration"],
                                       xj = data[, "Drug B Concentration"]),
                           start = start,
                           lower = lower,
                           upper = upper,
                           algorithm = "port",
                           control = nls.control(warnOnly = TRUE));
    fitResidual <- fit$m$deviance()},
        TRUE)
    
    minResidual <- min(c(startResidual, fitResidual))
    
    if (minResidual == startResidual)
    {
      katzirModelFit <- start
    }
    else
    {
      katzirModelFit <- fit$m$getPars()
    }
  }
  else
  {
    stop("modelType not yet implemented.")
  }
  
  katzirResidual <- minResidual
  
  return(list(fit = katzirModelFit, residual = katzirResidual))
}

analyzeInteraction <- function (drugs, cellLine, tissueType, dataset, modelType)
{
  maxLogN <- 2 # = log(n)
  minLogN <- -maxLogN
  load(paste0(tissueType, "_", dataset, "_twoDrugData.RData"))
  hillModelFits <- as.list(rep(NA, length(drugs)))
  names(hillModelFits) <- unlist(drugs)
  hillModelFitConfidenceIntervals <- hillModelFits
  
  if (modelType %in% c("bidirectional", "unidirectional"))
  {
    load(paste0(tissueType, "_", dataset, "_singleDrugData.RData"))
    singleDrugDoseResponseData <- hillModelFits
    
    for (drug in drugs)
    {
      singleDrugDoseResponseData[[drug]] <- extractDrugCellLineCombinationData(drugs = as.list(drug),
                                                                               cellLine = cellLine,
                                                                               allData = singleDrugData,
                                                                               dataType = "singleDrug")
      minLogD0 <- log(min(singleDrugDoseResponseData[[drug]][, "Concentration"])) - 1
      maxLogD0 <- log(max(singleDrugDoseResponseData[[drug]][, "Concentration"])) + 1
      hillModelFits[[drug]] <- nls(formula = y ~ Hill(x, logD0, logN),
                                   data = list(y = singleDrugDoseResponseData[[drug]][, "Viability"],
                                               x = singleDrugDoseResponseData[[drug]][, "Concentration"]),
                                   start = list(logD0 = log(singleDrugDoseResponseData[[drug]][which.min(abs(singleDrugDoseResponseData[[drug]][, "Viability"] - 0.5)), "Concentration"]),
                                                logN = 0),
                                   lower = list(logD0 = minLogD0, n = minLogN),
                                   upper = list(logD0 = maxLogD0, n = maxLogN),
                                   algorithm = "port",
                                   control = nls.control(warnOnly = TRUE))
      hillModelFitConfidenceIntervals[[drug]] <- matrix(c(-Inf, minLogN, Inf, maxLogN),
                                                        nrow = 2,
                                                        dimnames = list(c("logD0", "logN"), c("2.5%", "97.5%")))
      try(expr = {hillModelFitConfidenceIntervals[[drug]] <- confint(hillModelFits[[drug]])},
          TRUE)
      
      if (is.na(hillModelFitConfidenceIntervals[[drug]]["logD0", "2.5%"]))
      {
        hillModelFitConfidenceIntervals[[drug]]["logD0", "2.5%"] <- minLogD0
      }
      
      if (is.na(hillModelFitConfidenceIntervals[[drug]]["logD0", "97.5%"]))
      {
        hillModelFitConfidenceIntervals[[drug]]["logD0", "97.5%"] <- maxLogD0
      }
      
      if (is.na(hillModelFitConfidenceIntervals[[drug]]["logN", "2.5%"]))
      {
        hillModelFitConfidenceIntervals[[drug]]["logN", "2.5%"] <- minLogN
      }
      
      if (is.na(hillModelFitConfidenceIntervals[[drug]]["logN", "97.5%"]))
      {
        hillModelFitConfidenceIntervals[[drug]]["logN", "97.5%"] <- maxLogN
      }
    }
  }
  
  drugPairs <- list()
  
  for (i in 1:(length(drugs) - 1))
  {
    for (j in (i + 1):length(drugs))
    {
      drugPairs <- c(drugPairs, paste(drugs[[i]], "&", drugs[[j]]))
    }
  }
  
  katzirModelFits <- as.list(rep(NA, length(drugPairs)))
  names(katzirModelFits) <- unlist(drugPairs)
  katzirResiduals <- katzirModelFits
  twoDrugDoseResponseData <- katzirModelFits
  
  for (drugPair in drugPairs)
  {
    drugi <- unlist(strsplit(drugPair, " & "))[1]
    drugj <- unlist(strsplit(drugPair, " & "))[2]
    twoDrugDoseResponseData[[drugPair]] <- extractDrugCellLineCombinationData(drugs = as.list(c(drugi, drugj)),
                                                                              cellLine = cellLine,
                                                                              allData = twoDrugData,
                                                                              dataType = "twoDrug")
    result <- fitModelToData(data = twoDrugDoseResponseData[[drugPair]],
                             modelType = modelType,
                             hillModelFits = hillModelFits,
                             hillModelFitConfidenceIntervals = hillModelFitConfidenceIntervals)
    katzirModelFits[[drugPair]] <- result[["fit"]]
    katzirResiduals[[drugPair]] <- result[["residual"]]
  }
  
  analysis <- multiDrugCellLineCombinationAnalysis(drugs = as.list(drugs),
                                                   cellLine = cellLine,
                                                   tissueType = tissueType,
                                                   dataset = dataset,
                                                   hillModelFits = hillModelFits,
                                                   hillModelFitConfidenceIntervals = hillModelFitConfidenceIntervals,
                                                   katzirModelFits = katzirModelFits,
                                                   katzirResiduals = katzirResiduals)
  
  return(analysis)
}

refineKatzirStart <- function(xi, xj, y, start, lower, upper, missingParams = NULL)
{
  for (param in missingParams)
  {
    start[[param]] <- 0
    upper[[param]] <- 0
    lower[[param]] <- 0
  }
  
  start <- katzirMeshOptimization(xi,
                                  xj,
                                  y,
                                  unlist(start),
                                  unlist(lower),
                                  unlist(upper),
                                  c(0.5, 0.5, 0.5, 0.5, 0.25, 0.25),
                                  (names(start) %in% missingParams))
  start <- as.list(start)
  names(start) <- names(upper)
  
  if (length(missingParams) > 0)
  {
    start <- start[-which(names(start) %in% missingParams)]
  }
  
  return(start)
}

bootstrapParamCIs <- function(data,
                              modelType,
                              fittedParams,
                              hillModelFits,
                              hillModelFitConfidenceIntervals)
{
  nTrials <- 40
  simulatedData <- data
  
  if (!("logAij" %in% names(fittedParams)))
  {
    fittedParams[["logAij"]] <- 0
  }
  
  if (!("logAji" %in% names(fittedParams)))
  {
    fittedParams[["logAji"]] <- 0
  }
  
  residuals <- data[, "Viability"] - Katzir2(xi = data[, "Drug A Concentration"],
                                             xj = data[, "Drug B Concentration"],
                                             logD0i = fittedParams[["logD0i"]],
                                             logD0j = fittedParams[["logD0j"]],
                                             logNi = fittedParams[["logNi"]],
                                             logNj = fittedParams[["logNj"]],
                                             logAij = fittedParams[["logAij"]],
                                             logAji = fittedParams[["logAji"]])
  bootstrappedParams <- matrix(0, nrow = nTrials, ncol = length(fittedParams))
  colnames(bootstrappedParams) <- names(fittedParams)
  
  for (i in seq_len(nTrials))
  {
    simulatedData[, "Viability"] <- data[, "Viability"] + sample(x = residuals,
                                                                 size = nrow(data),
                                                                 replace = TRUE)
    simulatedParams <- fitModelToData(data = simulatedData,
                                      modelType = modelType,
                                      hillModelFits = hillModelFits,
                                      hillModelFitConfidenceIntervals = hillModelFitConfidenceIntervals)
    
    for (j in names(simulatedParams$fit))
    {
      bootstrappedParams[i, j] <- simulatedParams$fit[[j]]
    }
  }
  
  paramCIs <- matrix(NA, nrow = length(fittedParams), ncol = 2)
  rownames(paramCIs) <- names(fittedParams)
  colnames(paramCIs) <- c("2.5%", "97.5%")
  
  for (param in names(fittedParams))
  {
    sortedParamVals <- sort(bootstrappedParams[, param], decreasing = FALSE)
    paramCIs[param, "2.5%"] <- sortedParamVals[round(nTrials * 0.025) + 1]
    paramCIs[param, "97.5%"] <- sortedParamVals[round(nTrials * 0.975)]
  }
  
  return(paramCIs)
}