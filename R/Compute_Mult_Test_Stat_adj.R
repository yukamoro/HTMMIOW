#' Quantify indirect effect using inverse odds weighting adjusted for confounders
#'
#' This function performs mediation analysis using inverse odds weighting to
#' compute the indirect effect of multiple mediators with the presence of confounders.
#'
#' @param exposure vector of exposure variable
#' @param outcome vector of outcome variable
#' @param confounder matrix of covariates
#' @param mediator matrix of mediator variables
#' @return The indirect effect of multiple mediators
#' @export

Compute_Mult_Test_Stat_adj = function(exposure, outcome, confounder,mediator){
  #Fit exposure on mediators
  fit1 = glm(exposure ~ mediator + confounder,
             family = "binomial")
  #Compute inverse odds weights
  predprob = predict(fit1, type = "response")
  predprob = ifelse(predprob>1, 1, predprob)
  invodds = ((1-predprob)/predprob)
  invodd_weights = ifelse(exposure == 0, 1, invodds)
  #Fit outcome on exposure with weights for direct effect
  fit2 = glm(outcome ~ exposure + confounder,
             family = "quasipoisson",
             weights = invodd_weights)
  #Fit outcome on exposure for total effect
  fitte = glm(outcome ~ exposure + confounder,
              family = "quasipoisson")
  #Calculate indirect effect
  ie = fitte$coefficients[2] - fit2$coefficients[2]
  return(ie)
}
