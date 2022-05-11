Hill <- function(x, logD0, logN)
{
  return (1 / (1 + exp(exp(logN) * (log(x) - logD0))))
}

effectiveConcentrationFactor <- function(otherDrugEffectiveConcentrations, otherDrugLogD0, otherDrugLogA)
{
  D0 <- exp(otherDrugLogD0)
  A <- exp(otherDrugLogA)
  return((otherDrugEffectiveConcentrations + D0) / (A * otherDrugEffectiveConcentrations + D0))
}

effectiveConcentration <- function(myConcentrations, otherEffectiveConcentrations, otherLogD0s, otherLogAs)
{
  myEffectiveLogConcentrations <- log(myConcentrations)
  otherD0s <- exp(otherLogD0s)
  
  for (otherDrug in names(otherLogD0s))
  {
    myEffectiveLogConcentration <- myEffectiveLogConcentration - log(effectiveConcentrationFactor(otherEffectiveConcentrations[, otherDrug],
                                                                                                  otherDrugD0s[otherDrug],
                                                                                                  otherLogAs[otherDrug]))
  }
  
  return (exp(myEffectiveLogConcentration))
}

effectiveConcentration2 <- function(myConcentrations, otherConcentrations, myLogD0, otherLogD0, logAOtherOnMine, logAMineOnOther)
{
  myD0 <- exp(myLogD0)
  otherD0 <- exp(otherLogD0)
  AMineOnOther <- exp(logAMineOnOther)
  AOtherOnMine <- exp(logAOtherOnMine)
  A <- AMineOnOther * otherD0 + AOtherOnMine * otherConcentrations
  B <- myD0 * otherD0 + AOtherOnMine * otherConcentrations * myD0 -
    AMineOnOther * myConcentrations * otherD0 - myConcentrations * otherConcentrations
  C <- - myD0 * otherD0 * myConcentrations - myConcentrations * otherConcentrations * myD0
  
  return ((-B + sqrt(B ^ 2 - 4 * A * C)) / (2 * A))
}

effectiveConcentrationEstimateErrors <- function(concentrations, logD0s, logAs)
{
  errorEstimates <- rep(NA, length(concentrations))
  browser()
  for (i in seq_len(concentrations))
  {
    errorEstimates[i] <- concentrations[i] - effectiveConcentration(myConcentrations = concentrations[i],
                                                                    otherEffectiveConcentrations = concentrations[-i],
                                                                    otherLogD0s = logD0s[-i],
                                                                    otherLogAs = logAs[i, -i])
  }
  
  return(errorEstimates)
}

Katzir <- function(x, logD0s, logNs, logAs)
{
  if (length(logD0s) == 2)
  {
    return (Katzir2(xi = x[, 1],
                    xj = x[, 2],
                    logD0i = logD0s[1],
                    logD0j = logD0s[2],
                    logNi = logNs[1],
                    logNj = logNs[2],
                    logAij = logAs[1, 2],
                    logAji = logAs[2, 1]))
  }
  browser()
  effectiveConcentrations <- apply(x,
                                   1,
                                   function(row){return(nleqslv(x = row,
                                                                fn = effectiveConcentrationEstimateErrors(concentrations = row,
                                                                                                          logD0s = logD0s,
                                                                                                          logAs = logAs)))})
  viabilities <- apply(x,
                       1,
                       function(row){return(prod(sapply(row,
                                                        function(col){return(Hill(col,
                                                                                  logD0s[col],
                                                                                  logNs[col]))})))})
  
  return(viabilities)
}

Katzir2 <- function(xi, xj, logD0i, logD0j, logNi, logNj, logAij, logAji)
{
  return (Hill(effectiveConcentration2(xi,
                                       xj,
                                       logD0i,
                                       logD0j,
                                       logAij,
                                       logAji),
               logD0i,
               logNi) *
            Hill(effectiveConcentration2(xj,
                                         xi,
                                         logD0j,
                                         logD0i,
                                         logAji,
                                         logAij),
                 logD0j,
                 logNj))
}

KatzirUnidirectional <- function(x, logD0s, logNs, logAs)
{
  drugDependencyHierarchy <- order(apply(logAs,
                                         2,
                                         function(x){return(sum(x != 0))}),
                                   decreasing = FALSE)
  effectiveX <- x
  viability <- 1
  
  for (i in drugDependencyHierarchy)
  {
    for (j in seq_along(x))
    {
      if (logAs[j, i] != 0)
      {
        effectiveX[i] <- effectiveConcentration2(myConcentrations = x[i],
                                                 otherConcentrations = x[j],
                                                 myLogD0 = logD0s[i],
                                                 otherLogD0 = logD0s[j],
                                                 logAOtherOnMine = logAs[j, i],
                                                 logAMineOnOther = 0)
      }
    }
    
    viability <- viability * Hill(effectiveX[i], logD0s[i], logNs[i])
  }
  
  return(viability)
}