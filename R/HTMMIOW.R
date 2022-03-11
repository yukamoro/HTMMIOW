


library(dplyr)
library(ggplot2)
library(vegan)
library(compositions)
library(umap)
library(foreach)
library(doParallel)

registerDoParallel(20)
Compute_Mult_Test_Stat = function(exposure, outcome, mediator){
  #Fit exposure on outcome
  #mediator = meddata[,1]
  fit1 = glm(exposure ~ mediator,
             family = "binomial")
  #Compute an IORW/IOW weight
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
  #summary(fit2)
  #Calculate indirect effect (mediation)
  ie = fitte$coefficients[2] - fit2$coefficients[2]
  ie
}

Compute_Mult_Test_Stat_adj = function(exposure, outcome, confounder,mediator){
  #Fit exposure on outcome
  #mediator = meddata[,1]
  fit1 = glm(exposure ~ mediator + confounder,
             family = "binomial")
  #Compute an IORW/IOW weight
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
  #summary(fit2)
  #Calculate indirect effect (mediation)
  ie = fitte$coefficients[2] - fit2$coefficients[2]
  ie
}


#########################Perform HTMMIOW##########################

Perform_HTMMIOW = function(exposure, outcome, confounder = NULL, mediator, umapcomponents=2){
  #print("Tranforming Microbiome Data...")
  #Impute 0 values
  mediator_imp = mediator
  mediator_imp[mediator_imp==0] = 0.5
  #clr transform
  clr_trans = clr(mediator_imp)
  #PCA
  pca_mat  = prcomp(clr_trans)
  #UMAP
  #set.seed(100)
  umap_mat = umap(pca_mat$x, n_components = umapcomponents)
  transmedmat = umap_mat$layout

  #Hypothesis test
  #print("Performing Hypothesis Test...")
  if(is.null(confounder)){
    #Test Statistic
    tk = Compute_Mult_Test_Stat(exposure, outcome, transmedmat)
    #Compute Null
    B = 2000
    #set.seed(100)
    perm_null = replicate(B, Compute_Mult_Test_Stat(exposure, outcome, transmedmat[(sample(c(1:nrow(transmedmat)), nrow(transmedmat), replace = FALSE)),]))
    pk_2 = (sum(abs(perm_null) > abs(tk))/B)
  }
  else{
    #Test Statistic
    tk = Compute_Mult_Test_Stat_adj(exposure, outcome, confounder, transmedmat)
    #Compute Null
    B = 2000
    #set.seed(100)
    perm_null = replicate(B, Compute_Mult_Test_Stat_adj(exposure, outcome, confounder, transmedmat[(sample(c(1:nrow(transmedmat)), nrow(transmedmat), replace = FALSE)),]))
    pk_2 = (sum(abs(perm_null) > abs(tk))/B)
  }

  output = list(tk = tk,
                pk_2 = pk_2,
                perm_null = perm_null)
  return(output)
}

