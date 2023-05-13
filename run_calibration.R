source("support_function.R")
load("imputation_data.Rdata")


expit = function(x){exp(x)/(1 + exp(x))}
library(R.utils)
library(purrr)
library(doParallel)
library(foreach)
library(adaptMCMC)
args = commandArgs(trailingOnly=TRUE, asValues = TRUE, 
                   defaults = list(n_cores = 10,
                                n_chains = 100,
                                output_path = paste0("mcmc_runs/", Sys.Date(), "/"),
                                N.mcmc = 100,
                                N.samples = 10,
                                N.samples.outcome = 20,
                                batch.length = 10,
                                burnin = 100))

with(args,{
  if(!dir.exists(output_path)){
    dir.create(output_path)
  }
  N.reps = burnin + batch.length * N.mcmc
  samples.to.use = c(
    rep(0, 
        burnin), 
    rep(
      c(
        rep(0, 
            batch.length - 1),
        1), 
      N.mcmc)
  )
  

  sample.outcome = function(mcmc.samples, params.f, name, ...){
    n.params = nrow(mcmc.samples)
    N_iter = n.params * N.samples.outcome
    simulated.outcomes = matrix(nrow = N_iter, ncol = 6)
    
    foreach(i = seq(N_iter)) %do%{
      j = ceiling(i/N.samples.outcome)
      params = params.f(mcmc.samples[j,])
      simulated.outcomes[i,1:2] = FEMPREP.Likelihood(params, 1, return.simulation = T,...)
      simulated.outcomes[i,3:4] = VOICE.Likelihood(params, 1, return.simulation = T,...)
      simulated.outcomes[i,5:6] = PARTNERSPREP.Likelihood(params, 1, return.simulation = T, ...)
    }
    
    colnames(simulated.outcomes) <- c("FEMPREP_POS", "FEMPREP_NEG", 
                                   "VOICE_POS", "VOICE_NEG",
                                   "PARTNERSPREP_POS", "PARTNERSPREP_NEG")
    data.table::data.table(simulated.outcomes)
    #save(simulated.outcomes, file = paste0(output_path, name, "_outcomes.Rdata"))
  }
  
  parallel.mcmc = function(LL, init, init.scale, name, params.f, ...){
    out.MCMC = foreach(i = seq(n_chains))%dopar%{
      .out = MCMC(LL, N.reps, init, scale = init.scale, K=N.samples, showProgressBar = T)
      
      .out$outcomes = sample.outcome(.out$samples, params.f, name, ...)
      
      .out
    }
    
    save(out.MCMC, file = paste0(output_path, name, ".Rdata"))
  }
  
  registerDoParallel(cores = n_cores)
  
  
  parallel.mcmc(Log.Likelihood.PLASMA, c(0, 0), c(0.1, 0.1), "plasma", params.f.plasma, plasma.only = T)
  parallel.mcmc(Log.Likelihood.PLASMA2, 0, 0.1, "plasma2", params.f.plasma2, plasma.only = T)
  parallel.mcmc(Log.Likelihood.DBS, 650, 100, "DBS",  params.f.DBS, f = hazard.ratio.by.dbs)
  parallel.mcmc(Log.Likelihood.DBS2, c(650,0), c(100, 0.1), "DBS2", params.f.DBS2,  f = hazard.ratio.by.dbs2)
  parallel.mcmc(Log.Likelihood.DBS3, c(650, -3), c(100,0.1), "DBS3",  params.f.DBS3, f = hazard.ratio.by.dbs3)
})