
#' Perform HT-MMIOW
#'
#' This function performs HT-MMIOW on user-provided exposure, outcome,
#' confounder, and mediator data.
#'
#'
#' @param exposure vector of exposure variable
#' @param outcome vector of outcome variable
#' @param confounder matrix of relavant covariates
#' @param mediator matrix of integers that represent abundance
#' @param umapcomponents number of components for UMAP to reduce to
#' @param B number of replicates for hypothesis test
#' @return A list of observed test statistic, p-value, and test statistics under the null hypothesis
#' @export
#' @importFrom compositions clr
#' @importFrom umap umap
#' @importFrom stats prcomp
#' @examples
#' Perform_HTMMIOW(exposure = sim_data$sim_exposure,
#' outcome = sim_data$sim_outcome,
#' mediator = sim_data$sim_microbiome,
#' umapcomponents = 2)


Perform_HTMMIOW = function(exposure,
                           outcome,
                           confounder = NULL,
                           mediator,
                           umapcomponents=2,
                           B = 2000){
  #print("Tranforming Microbiome Data...")
  #Impute 0 values with pseudocount
  mediator_imp = mediator
  mediator_imp[mediator_imp==0] = 0.5
  #clr transform
  clr_trans = clr(mediator_imp)
  #SVD
  pca_mat  = prcomp(clr_trans)
  #UMAP
  umap_mat = umap(pca_mat$x, n_components = umapcomponents)
  transmedmat = umap_mat$layout

  #Perform hypothesis test
  if(is.null(confounder)){
    #Observed Test Statistic (Indirect Effect)
    tk = Compute_Mult_Test_Stat(exposure, outcome, transmedmat)
    #Compute statistics under the null hypothesis
    perm_null = replicate(B,
                          Compute_Mult_Test_Stat(exposure,
                                                 outcome,
                                                 transmedmat[(sample(c(1:nrow(transmedmat)), nrow(transmedmat), replace = FALSE)),]))
    #Calculate P-value
    pk_2 = (sum(abs(perm_null) > abs(tk))/B)
  }
  else{
    #Observed Test Statistic (Indirect Effect)
    tk = Compute_Mult_Test_Stat_adj(exposure, outcome, confounder, transmedmat)
    #Compute statistics under the null hypothesis
    perm_null = replicate(B,
                          Compute_Mult_Test_Stat_adj(exposure,
                                                     outcome,
                                                     confounder,
                                                     transmedmat[(sample(c(1:nrow(transmedmat)), nrow(transmedmat), replace = FALSE)),]))
    #Calculate P-value
    pk_2 = (sum(abs(perm_null) > abs(tk))/B)
  }
  #Create Output
  output = list(tk = tk,
                pk_2 = pk_2,
                perm_null = perm_null)
  return(output)
}

