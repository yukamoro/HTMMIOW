#Simulate data for HT-MMIOW

library(SparseDOSSA2)

Simulate_data = function(n_samples,
                         n_microbes,
                         ef_mede,
                         perc_spike,
                         es_exposure,
                         es_mediators,
                         num_meds){


  checkMedMicrobe = 0
  while(checkMedMicrobe==0){
    Stool_simulation <- SparseDOSSA2(template = "Stool",  # choose from "Stool", "Vaginal" or "IBD"
                                     new_features = TRUE,  # should new features be simulated
                                     n_sample = n_samples,  # number of samples to simulate
                                     n_feature = n_microbes,  # number of features to simulate (when 'new_features = TRUE')
                                     verbose = FALSE, # return detailed info
                                     spike_metadata = "abundance", #association with metadata will be spiked in
                                     metadata_effect_size = ef_mede, #effect size of metadata
                                     perc_feature_spiked_metadata = perc_spike)  #Percent of features associated

    #Saving output to variables
    sim_microbiome = t(Stool_simulation$simulated_data)
    sim_microbiome_relabun = t(Stool_simulation$simulated_matrices$rel)
    sim_exposure = Stool_simulation$spike_metadata$metadata_matrix[,2]
    #Microbiome mediators
    spiked_microbes = Stool_simulation$spike_metadata$feature_metadata_spike_df
    med_microbes = spiked_microbes$feature_spiked[which(spiked_microbes$metadata_datum==2)]

    #Remove spiked microbes that are majority zeros
    med_microbes = med_microbes[which(colSums(sim_microbiome[,med_microbes])>0)]

    #Check if # of selected microbes is larger than number of mediators of interest
    checkMedMicrobe = ifelse(length(med_microbes)>=num_meds, 1, 0)
  }
  #print("E and M Simulation Done")
  med_microbes = med_microbes[sample(length(med_microbes), num_meds, replace = FALSE)]
  #Create vector of outcome~mediator coefficients
  es_mediators_vec = rep(es_mediators,length(med_microbes))

  #print(med_microbes)
  #colSums(sim_microbiome[,med_microbes])

  #Simulate outcome
  sim_microbiome_relabun_2 = sim_microbiome_relabun
  sim_microbiome_relabun_2[sim_microbiome_relabun_2==0] = 5e-20
  #Predict outcome using parameters
  estimate = es_exposure*sim_exposure + scale(log(sim_microbiome_relabun_2[,med_microbes]))%*%es_mediators_vec + rnorm(n_samples, 0, 1)
  sim_outcome = ifelse(exp(estimate)/(1+exp(estimate)) >0.5, 1, 0)
  #print("Outcome Simulation Done")
  #table(sim_outcome, sim_exposure)

  #ht = Perform_HTMMIOW(exposure = sim_exposure, outcome = sim_outcome, mediator = sim_microbiome)
  #ht$pk

  #Calculate R2
  x = scale(log(sim_microbiome_relabun_2[,med_microbes]))[,1]
  R2 = ((n_samples*sum(x * estimate) - sum(x)*sum(estimate)) / (sqrt(n_samples*sum(x^2)-(sum(x)^2)) * sqrt(n_samples*sum(estimate^2)-(sum(estimate)^2))))^2
  R2
  #coef = 0.5: R2 for each taxa = 0.05479735
  #coef = 1: R2 for each taxa = 0.3020513
  #coef = 1: R2 for each taxa = 0.8800833



  output = list(sim_microbiome = sim_microbiome,
                sim_exposure = sim_exposure,
                sim_outcome = sim_outcome)
  return(output)
}


set.seed(10)
#Simulate microbiome and exposure using SparseDOSSA2 from BioBakery
sim_data = Simulate_data(n_samples = 300,
                         n_microbes = 400,
                         ef_mede = 5,
                         perc_spike = 0.5,
                         es_exposure = 3,
                         es_mediators = 5,
                         num_meds = 10)
