multiDrugCellLineCombinationAnalysis_lowData <- setRefClass("multiDrugCellLineCombinationAnalysis",
                                                            fields = list(drugs = "list",
                                                                          cellLine = "character",
                                                                          twoDrugDoseResponseData = "list",
                                                                          katzirModelFits = "list"),
                                                            methods = list(
                                                              drugPairs = function()
                                                              {
                                                                drugPairs <- list()
                                                                
                                                                for (i in 1:(length(unlist(drugs)) - 1))
                                                                {
                                                                  for (j in (i + 1):length(unlist(drugs)))
                                                                  {
                                                                    drugPairs <- c(drugPairs, paste(drugs[[1]][[i]], "&", drugs[[1]][[j]]))
                                                                  }
                                                                }
                                                                
                                                                return(drugPairs)
                                                              },
                                                              plotResults = function()
                                                              {
                                                                for (drugPair in drugPairs())
                                                                {
                                                                  xs <- sort(unique(twoDrugDoseResponseData[[drugPair]][, "Drug i Concentration"]))
                                                                  ys <- sort(unique(twoDrugDoseResponseData[[drugPair]][, "Drug j Concentration"]))
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
                                                                        xlab = "Drug i Log Concentration (uM)",
                                                                        ylab = "Drug j Log Concentration (uM)",
                                                                        col = colourScale)
                                                                  points(x = log(twoDrugDoseResponseData[[drugPair]][, "Drug i Concentration"]),
                                                                         y = log(twoDrugDoseResponseData[[drugPair]][, "Drug j Concentration"]),
                                                                         pch = 21,
                                                                         col = "black",
                                                                         bg = colourScale[ceiling(twoDrugDoseResponseData[[drugPair]][, "Viability"] *
                                                                                                    length(colourScale))])
                                                                }
                                                              }))

analyzeInteraction_lowData <- function (drugs, cellLine, twoDrugData)
{
  maxLogD0 <- 18 # = log concentration of water in water in umol
  maxLogN <- 2 # = log(n)
  minLogD0 <- -41 # = log concentration of 1 atom in 1 L water in umol
  minLogN <- -maxLogN
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
  maxLogA <- 5 # = log(1 + a)
  minLogA <- -maxLogA
  twoDrugConcentrationColumns <- c("Drug A Concentration", "Drug B Concentration")
  twoDrugDoseResponseData <- katzirModelFits
  twoDrugViabilityColumns <- grep("Viability", colnames(twoDrugData))
  
  for (drugPair in drugPairs)
  {
    drugi <- unlist(strsplit(drugPair, " & "))[1]
    drugj <- unlist(strsplit(drugPair, " & "))[2]
    combinationRowsAB <- which(twoDrugData[, "Drug A"] == drugi &
                                 twoDrugData[, "Drug B"] == drugj &
                                 twoDrugData[, "Cell Line"] == cellLine)
    combinationRowsBA <- which(twoDrugData[, "Drug A"] == drugj &
                                 twoDrugData[, "Drug B"] == drugi &
                                 twoDrugData[, "Cell Line"] == cellLine)
    twoDrugDoseResponseData[[drugPair]] <- rbind(twoDrugData[combinationRowsAB, ],
                                                 twoDrugData[combinationRowsBA, c("Cell Line",
                                                                                  "Drug B",
                                                                                  "Drug A",
                                                                                  "Drug B Concentration",
                                                                                  "Drug A Concentration",
                                                                                  "Viability")])
    colnames(twoDrugDoseResponseData[[drugPair]]) <- c("Cell Line",
                                                       "Drug A",
                                                       "Drug B",
                                                       "Drug A Concentration",
                                                       "Drug B Concentration",
                                                       "Viability")
    start <- list(logD0i = log(twoDrugDoseResponseData[[drugPair]][which.min(abs(twoDrugDoseResponseData[[drugPair]][, "Viability"] - 0.25)), "Drug A Concentration"]),
                  logD0j = log(twoDrugDoseResponseData[[drugPair]][which.min(abs(twoDrugDoseResponseData[[drugPair]][, "Viability"] - 0.25)), "Drug B Concentration"]),
                  logAji = 0)
    lower <- list(logD0i = minLogD0,
                  logD0j = minLogD0,
                  logAji = minLogA)
    upper <- list(logD0i = maxLogD0,
                  logD0j = maxLogD0,
                  logAji = maxLogA)
    start <- refineKatzirStart(y = twoDrugDoseResponseData[[drugPair]][, "Viability"],
                               xi = twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"],
                               xj = twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"],
                               start = start,
                               lower = lower,
                               upper = upper,
                               missingParams = c("logNi", "logNj", "logAij"))
    fit1 <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, 0, 0, 0, logAji),
                data = list(y = twoDrugDoseResponseData[[drugPair]][, "Viability"],
                            xi = twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"],
                            xj = twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"]),
                start = start,
                lower = lower,
                upper = upper,
                algorithm = "port",
                control = nls.control(warnOnly = TRUE))
    
    start <- list(logD0i = log(median(twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"])),
                  logD0j = log(median(twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"])),
                  logAij = 0)
    lower <- list(logD0i = minLogD0,
                  logD0j = minLogD0,
                  logAij = minLogA)
    upper <- list(logD0i = maxLogD0,
                  logD0j = maxLogD0,
                  logAij = maxLogA)
    start <- refineKatzirStart(y = twoDrugDoseResponseData[[drugPair]][, "Viability"],
                               xi = twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"],
                               xj = twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"],
                               start = start,
                               lower = lower,
                               upper = upper,
                               missingParams = c("logNi", "logNj", "logAji"))
    startResidual <- sum((y - Katzir2(xi,
                                      xj,
                                      start[["logD0i"]],
                                      start[["logD0j"]],
                                      start[["logNi"]],
                                      start[["logNj"]],
                                      start[["logAij"]],
                                      start[["logAji"]])) ^ 2)
    fit2 <- nls(formula = y ~ Katzir2(xi, xj, logD0i, logD0j, 0, 0, logAij, 0),
                data = list(y = twoDrugDoseResponseData[[drugPair]][, "Viability"],
                            xi = twoDrugDoseResponseData[[drugPair]][, "Drug A Concentration"],
                            xj = twoDrugDoseResponseData[[drugPair]][, "Drug B Concentration"]),
                start = start,
                lower = lower,
                upper = upper,
                algorithm = "port",
                control = nls.control(warnOnly = TRUE))
    
    if (fit1$m$deviance() <= fit2$m$deviance())
    {
      katzirModelFits[[drugPair]] <- fit1
    }
    else
    {
      katzirModelFits[[drugPair]] <- fit2
    }
    
    if (fit$m$deviance() < startResidual)
    {
      katzirModelFits[[drugPair]] <- fit$m$getPars()
    }
    else
    {
      katzirModelFits[[drugPair]] <- start
    }
  }
  
  analysis <- multiDrugCellLineCombinationAnalysis_lowData(drugs = list(drugs),
                                                           cellLine = cellLine,
                                                           twoDrugDoseResponseData = twoDrugDoseResponseData,
                                                           katzirModelFits = katzirModelFits)
  print(katzirModelFits[[drugPair]]$m$deviance())
  return(analysis)
}