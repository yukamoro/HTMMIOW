#' Quantify indirect effect using inverse odds weighting
#'
#' This function performs mediation analysis using inverse
#' odds weighting to compute the indirect effect of multiple mediators.
#'
#' @param exposure vector of exposure variable
#' @param outcome vector of outcome variable
#' @param mediator matrix of mediator variables
#' @return The indirect effect of multiple mediators
#' @export

Compute_Mult_Test_Stat = function(exposure, outcome, mediator){
  #Fit exposure on mediators
  fit1 = glm(exposure ~ mediator,
             family = "binomial")
  #Compute inverse odds weights
  predprob = predict(fit1, type = "response")
  predprob = ifelse(predprob>1, 1, predprob)
  invodds = ((1-predprob)/predprob)
  invodd_weights = ifelse(exposure == 0, 1, invodds)
  #Fit outcome on exposure with weights for direct effect
  fit2 = glm(outcome ~ exposure ,
             family = "quasipoisson",
             weights = invodd_weights)
  #Fit outcome on exposure for total effect
  fitte = glm(outcome ~ exposure ,
              family = "quasipoisson")
  #Calculate indirect effect
  ie = fitte$coefficients[2] - fit2$coefficients[2]
  return(ie)
}
