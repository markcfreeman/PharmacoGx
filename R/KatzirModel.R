.effectiveConcentrationFactor <- function(otherDrugEffectiveConcentrations, otherDrugLogD0, otherDrugLogA)
{
  D0 <- exp(otherDrugLogD0)
  A <- exp(otherDrugLogA)
  return((otherDrugEffectiveConcentrations + D0) / (A * otherDrugEffectiveConcentrations + D0))
}

.effectiveConcentration <- function(myConcentrations, otherEffectiveConcentrations, otherLogD0s, otherLogAs)
{
  myEffectiveLogConcentrations <- log(myConcentrations)
  otherD0s <- exp(otherLogD0s)
  
  for (otherDrug in names(otherLogD0s))
  {
    myEffectiveLogConcentration <- myEffectiveLogConcentration - log(PharmacoGx:::.effectiveConcentrationFactor(otherEffectiveConcentrations[, otherDrug],
                                                                                                                otherDrugD0s[otherDrug],
                                                                                                                otherLogAs[otherDrug]))
  }
  
  return (exp(myEffectiveLogConcentration))
}

.effectiveConcentration2 <- function(myConcentrations, otherConcentrations, myLogD0, otherLogD0, logAOtherOnMine, logAMineOnOther)
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

.effectiveConcentrationEstimateErrors <- function(concentrations, logD0s, logAs)
{
  errorEstimates <- rep(NA, length(concentrations))
  browser()
  for (i in seq_len(concentrations))
  {
    errorEstimates[i] <- concentrations[i] - PharmacoGx:::.effectiveConcentration(myConcentrations = concentrations[i],
                                                                                  otherEffectiveConcentrations = concentrations[-i],
                                                                                  otherLogD0s = logD0s[-i],
                                                                                  otherLogAs = logAs[i, -i])
  }
  
  return(errorEstimates)
}

#' Fits curves of the form * 
#'
#' @examples
#' x <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' conc <- matrix(rep(x, 3), ncol = 3)
#' EC50 <- c(1.5, -0.2, 0.3)
#' HS <- c(0.9, 1, 0.95)
#' A <- matrix(c(1, 0.98, 1.22, 0.45, 1, 1.32, 0.92, 1.08, 1), ncol = 3)
#' Zimmer(conc, EC50, HS, A)
#' 
#' @param conc `numeric` is a *
#' @param EC50 `numeric` is a *   
#' @param HS `numeric` is a *
#' @param A `numeric` is a *
#'
#' @export
#' 
#' @importFrom nleqslv nleqslv
Katzir <- function(x, logD0s, logNs, logAs)
{
  if (length(logD0s) == 2)
  {
    return (PharmacoGx:::.Katzir2(xi = x[, 1],
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
                                                                fn = PharmacoGx:::.effectiveConcentrationEstimateErrors(concentrations = row,
                                                                                                                        logD0s = logD0s,
                                                                                                                        logAs = logAs)))})
  viabilities <- apply(x,
                       1,
                       function(row){return(prod(sapply(row,
                                                        function(col){return(PharmacoGx:::.Hill(col,
                                                                                                c(exp(logNs[col]),
                                                                                                  0,
                                                                                                  exp(logD0s[col]))))})))})
  
  return(viabilities)
}

.Katzir2 <- function(xi, xj, logD0i, logD0j, logNi, logNj, logAij, logAji)
{
  return (PharmacoGx:::.Hill(PharmacoGx:::.effectiveConcentration2(xi,
                                                                   xj,
                                                                   logD0i,
                                                                   logD0j,
                                                                   logAij,
                                                                   logAji),
                             c(exp(logNi),
                               0,
                               exp(logD0i))) *
            PharmacoGx:::.Hill(PharmacoGx:::.effectiveConcentration2(xj,
                                                                     xi,
                                                                     logD0j,
                                                                     logD0i,
                                                                     logAji,
                                                                     logAij),
                               c(exp(logNj),
                                 0,
                                 exp(logD0j))))
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
        effectiveX[i] <- PharmacoGx:::.effectiveConcentration2(myConcentrations = x[i],
                                                               otherConcentrations = x[j],
                                                               myLogD0 = logD0s[i],
                                                               otherLogD0 = logD0s[j],
                                                               logAOtherOnMine = logAs[j, i],
                                                               logAMineOnOther = 0)
      }
    }
    
    viability <- viability * PharmacoGx:::.Hill(effectiveX[i],
                                                c(exp(logNs[i]),
                                                  0,
                                                  exp(logD0s[i])))
  }
  
  return(viability)
}