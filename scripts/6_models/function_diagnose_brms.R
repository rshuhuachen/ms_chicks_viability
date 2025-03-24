
diagnose <- function(fit, modelname){
  #get posteriors
  posterior <- as.array(fit)
  log_ps <- log_posterior(fit)
  nuts <- nuts_params(fit) #divergence
  #get only beta and sd
  betas <- variables(fit)[grep("b_", variables(fit))]
  sd <- variables(fit)[grep("sd_", variables(fit))]
  
  #global patterns in divergence
  diverge_beta <- mcmc_parcoord(posterior, np = nuts, pars= betas)
  
  #identify collinearity between parameters
  collin_beta <- mcmc_pairs(posterior, np = nuts, pars= betas)
  
  #traceplot
  trace_beta <- mcmc_trace(posterior, pars = betas, np = nuts)
  trace_sd <- mcmc_trace(posterior, pars = sd, np = nuts)
  
  #rhat
  rhat <- mcmc_rhat(brms::rhat(fit))
  
  #effective sample size
  neff <- mcmc_neff(neff_ratio(fit))
  
  #autocorrelation
  autocor_beta <- mcmc_acf(posterior, pars = betas)
  autocor_sd <- mcmc_acf(posterior, pars=sd)
  
  #quick glance results
  areas <- mcmc_areas(fit, pars=betas)
  
  #combine in list
  diagnosis <- list(diverge_beta = diverge_beta, 
                    collin_beta = collin_beta, 
                    trace_beta = trace_beta, 
                    trace_sd = trace_sd, 
                    rhat = rhat, 
                    neff = neff, 
                    autocor_beta = autocor_beta, 
                    autocor_sd = autocor_sd,
                    areas = areas)
  
  
  # save output
  pdf(file=paste0("output/models/diagnosis/", modelname, ".pdf"))
  print(diverge_beta)
  print(collin_beta)
  print(trace_beta)
  print(trace_sd)
  print(rhat)
  print(neff)
  print(autocor_beta)
  print(autocor_sd)
  print(areas)
  dev.off()
  
  # add to summary
  return(diagnosis)
}
