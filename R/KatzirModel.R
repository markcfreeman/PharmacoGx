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

#' Predicts the viability of cells exposed to combinations of different drugs at
#' different concentrations acting in accordance with the drug synergy model of
#' Zimmer et al. in their 2016 paper "Prediction of multidimensional drug dose
#' responses based on measurements of drug pairs".
#'
#' @examples
#' x <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' conc <- matrix(rep(x, 3), ncol = 3)
#' EC50 <- c(1.5, -0.2, 0.3)
#' HS <- c(0.9, 1, 0.95)
#' A <- matrix(c(1, 0.98, 1.22, 0.45, 1, 1.32, 0.92, 1.08, 1), ncol = 3)
#' Zimmer(conc, EC50, HS, A)
#' 
#' @param conc `numeric` is an m-by-n matrix where each row represents a
#' combination of drug concentrations acting simultaneously on a cell line.
#' @param EC50 `numeric` is a vector of length 'd' such that each element of the
#' vector is the EC50 of the drug whose concentrations are listed as the
#' corresponding element of 'conc'. Values in the vector correspond to the
#' 'D_0i' #' values of the Zimmer et al. model.
#' @param HS `numeric` is a vector of length 'd' such that each element of the
#' vector is the Hill slope of the drug whose concentrations are listed as the
#' corresponding element of 'conc'. Values in the vector correspond to the
#' 'n_i' values of the Zimmer et al. model.
#' @param a `numeric` is a d-by-d matrix such that the (i, j)th element of the
#' matrix is the 'a_ij' value of the corresponding drug pair in the Zimmer et
#' al. model. Note that entries along the main diagonal should equal 1.
#' @param unidirectional `logical` is TRUE if the simplifying assumption
#' described in Zimmer et al. wherein drug synergy is unidirectional should be
#' applied, and therefore that a_ji != 1 implies a_ij = 1.
#'
#' @export
#' 
#' @importFrom nleqslv nleqslv
Katzir <- function(conc, EC50, HS, a, unidirectional = FALSE)
{
  if (unidirectional)
  {
    if (length(logD0s) == 2)
    {
      return (PharmacoGx:::.Katzir2(xi = conc[[1]],
                                    xj = conc[[2]],
                                    logD0i = log(EC50[1]),
                                    logD0j = log(EC50[2]),
                                    logNi = log(HS[1]),
                                    logNj = log(HS[2]),
                                    logAij = log(a[1, 2]),
                                    logAji = log(a[2, 1])))
    }
    
    effectiveConcentrations <- apply(conc,
                                     1,
                                     function(x){return(nleqslv(x,
                                                                fn = PharmacoGx:::.effectiveConcentrationEstimateErrors(concentrations = concs,
                                                                                                                        logD0s = log(EC50),
                                                                                                                        logAs = log(a))))})
    viabilities <- lapply(1:length(conc),
                          function(x){return(prod(PharmacoGx:::.Hill(conc[[x]],
                                                                     c(HS[x],
                                                                       0,
                                                                       EC50[x]))))})
  }
  else
  {
    drugDependencyHierarchy <- order(apply(log(a),
                                           2,
                                           function(x){return(sum(x != 0))}),
                                     decreasing = FALSE)
    effectiveConc <- conc
    viability <- 1
    
    for (i in drugDependencyHierarchy)
    {
      for (j in seq_along(conc))
      {
        if (a[j, i] != 1)
        {
          effectiveConc[[i]] <- PharmacoGx:::.effectiveConcentration2(myConcentrations = conc[[i]],
                                                                      otherConcentrations = conc[[j]],
                                                                      myLogD0 = log(EC50[i]),
                                                                      otherLogD0 = log(EC50[j]),
                                                                      logAOtherOnMine = log(a[j, i]),
                                                                      logAMineOnOther = 0)
        }
      }
      
      viability <- viability * PharmacoGx:::.Hill(effectiveConc[i],
                                                  c(HS[i],
                                                    0,
                                                    EC50[i]))
    }
  }
  
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

#' Predicts the viability of cells exposed to combinations of different drugs at
#' different concentrations acting in accordance with the drug synergy model of
#' Zimmer et al. in their 2016 paper "Prediction of multidimensional drug dose
#' responses based on measurements of drug pairs", under 
#'
#' @examples
#' x <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' conc <- matrix(rep(x, 3), ncol = 3)
#' EC50 <- c(1.5, -0.2, 0.3)
#' HS <- c(0.9, 1, 0.95)
#' A <- matrix(c(1, 0.98, 1.22, 0.45, 1, 1.32, 0.92, 1.08, 1), ncol = 3)
#' ZimmerUnidirectional(conc, EC50, HS, A)
#' 
#' @param conc `numeric` is list of length 'd' such that each element of the
#' list is a vectors of drug of a different drug
#' @param EC50 `numeric` is a vector of length 'd' such that each element of the
#' vector is the EC50 of the drug whose concentrations are listed as the
#' corresponding element of 'conc'. Values in the vector correspond to the
#' 'D_0i' #' values of the Zimmer et al. model.
#' @param HS `numeric` is a vector of length 'd' such that each element of the
#' vector is the Hill slope of the drug whose concentrations are listed as the
#' corresponding element of 'conc'. Values in the vector correspond to the
#' 'n_i' values of the Zimmer et al. model.
#' @param a `matrix` is a d-by-d matrix such that the (i, j)th element of the
#' matrix is the 'a_ij' value of the corresponding drug pair in the Zimmer et
#' al. model. Note that entries along the main diagonal should equal 1.
#'
#' @export
#' 
#' @importFrom nleqslv nleqslv
KatzirUnidirectional <- function(conc, EC50, HS, a)
#KatzirUnidirectional <- function(x, logD0s, logNs, logAs)
{
  
  
  return(viability)
}