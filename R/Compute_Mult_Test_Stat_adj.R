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

  if(is.null(confounder)){
    #mediator = meddata[,1]
    fit1 = glm(exposure ~ mediator,
               family = "binomial")
    #Compute an IORW/IOW weight
    predprob = predict(fit1, type = "response")
    predprob = ifelse(predprob>1, 1, predprob)
    invodds = ((1-predprob)/predprob)
    invodd_weights = ifelse(exposure == 0, 1, invodds)

    #Check if dichotomous or continuous
    if(length(table(outcome))==2){
      #Fit outcome on exposure with weights for direct effect
      fit2 = glm(outcome ~ exposure ,
                 family = "quasipoisson",
                 weights = invodd_weights)
      #Fit outcome on exposure for total effect
      fitte = glm(outcome ~ exposure ,
                  family = "quasipoisson")
    } else{
      #Fit outcome on exposure with weights for direct effect
      fit2 = glm(outcome ~ exposure ,
                 weights = invodd_weights)
      #Fit outcome on exposure for total effect
      fitte = glm(outcome ~ exposure)
    }

  }
  else{
    #Fit exposure on outcome
    #mediator = meddata[,1]
    fit1 = glm(exposure ~ mediator + confounder,
               family = "binomial")
    #Compute an IORW/IOW weight
    predprob = predict(fit1, type = "response")
    predprob = ifelse(predprob>1, 1, predprob)
    invodds = ((1-predprob)/predprob)
    invodd_weights = ifelse(exposure == 0, 1, invodds)

    #Check if dichotomous or continuous
    if(length(table(outcome))==2){
      #Fit outcome on exposure with weights for direct effect
      fit2 = glm(outcome ~ exposure + confounder,
                 family = "quasipoisson",
                 weights = invodd_weights)
      #Fit outcome on exposure for total effect
      fitte = glm(outcome ~ exposure + confounder,
                  family = "quasipoisson")
    } else{
      #Fit outcome on exposure with weights for direct effect
      fit2 = glm(outcome ~ exposure + confounder,
                 weights = invodd_weights)
      #Fit outcome on exposure for total effect
      fitte = glm(outcome ~ exposure + confounder)
    }

  }


  #Calculate indirect effect (mediation)
  ie = fitte$coefficients[2] - fit2$coefficients[2]
  return(ie)
}
