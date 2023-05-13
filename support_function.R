# Model Definitions ####

hazard.ratio.by.dbs = function(dbs, params){
  with(as.list(params),{
    exp(-dbs/k.dbs)
  })
}

hazard.ratio.by.dbs2 = function(dbs, params){
  with(as.list(params),{
    hr.max * exp(-dbs/k.dbs)
  })
}

hazard.ratio.by.dbs3 = function(dbs, params){
  with(as.list(params),{
    hr.min + (1 - hr.min) * exp(-dbs/k.dbs)
  })
}

#Conversion of MCMC parameters to model parameters
params.f.DBS = function(x){
  list(k.dbs = x)
}

params.f.DBS2 = function(x){
  list(k.dbs = x[1], hr.max = exp(x[2]))
}

params.f.DBS3 = function(x){
  list(k.dbs = x[1], hr.min = expit(x[2]))
}

params.f.plasma = function(x){
  list(positive.hazard = expit(-x[1]), negative.hazard = exp(x[2]))
}

params.f.plasma2 = function(x){
  list(positive.hazard = expit(-x[1]), negative.hazard = 1)
}

# Likelihood Calculation ####

#K = number of simulations with each parameter set
Log.Likelihood.DBS = function(x, K=1, return.max = F){
  LL.VOICE = VOICE.Likelihood(list(k.dbs = x), K, return.max = return.max, f = hazard.ratio.by.dbs)
  LL.FEMPREP = FEMPREP.Likelihood(list(k.dbs = x), K, return.max = return.max, f = hazard.ratio.by.dbs)
  LL.PARTNERSPREP = PARTNERSPREP.Likelihood(list(k.dbs = x), K, return.max = return.max, f = hazard.ratio.by.dbs)
  LL.VOICE + LL.FEMPREP + LL.PARTNERSPREP
}

Log.Likelihood.DBS2 = function(x, K=1, return.max = F){
  LL.VOICE = VOICE.Likelihood(list(k.dbs = x[1], hr.max = exp(x[2])), K, return.max = return.max, f = hazard.ratio.by.dbs2)
  LL.FEMPREP = FEMPREP.Likelihood(list(k.dbs = x[1], hr.max = exp(x[2])), K, return.max = return.max, f = hazard.ratio.by.dbs2)
  LL.PARTNERSPREP = PARTNERSPREP.Likelihood(list(k.dbs = x[1], hr.max = exp(x[2])), K, return.max = return.max, f = hazard.ratio.by.dbs2)
  LL.VOICE + LL.FEMPREP + LL.PARTNERSPREP
}

Log.Likelihood.DBS3 = function(x, K=1, return.max = F){
  LL.VOICE = VOICE.Likelihood(list(k.dbs = x[1], hr.min = expit(x[2])), K, return.max = return.max, f = hazard.ratio.by.dbs3)
  LL.FEMPREP = FEMPREP.Likelihood(list(k.dbs = x[1], hr.min = expit(x[2])), K, return.max = return.max, f = hazard.ratio.by.dbs3)
  LL.PARTNERSPREP = PARTNERSPREP.Likelihood(list(k.dbs = x[1], hr.min = expit(x[2])), K, return.max = return.max, f = hazard.ratio.by.dbs3)
  LL.VOICE + LL.FEMPREP + LL.PARTNERSPREP
}

Log.Likelihood.PLASMA = function(x, K=1, return.max = F){
  params = list(positive.hazard = exp(x[1]), negative.hazard = exp(x[2]))
  LL.VOICE = VOICE.Likelihood(params, K, return.max = return.max, plasma.only = T)
  LL.FEMPREP = FEMPREP.Likelihood(params, K, return.max = return.max, plasma.only = T)
  LL.PARTNERSPREP = PARTNERSPREP.Likelihood(params, K, return.max = return.max, plasma.only = T)
  LL.VOICE + LL.FEMPREP + LL.PARTNERSPREP
}

Log.Likelihood.PLASMA2 = function(x, K=1, return.max = F){
  params = list(positive.hazard = exp(x[1]), negative.hazard = 1)
  LL.VOICE = VOICE.Likelihood(params, K, return.max = return.max, plasma.only = T)
  LL.FEMPREP = FEMPREP.Likelihood(params, K, return.max = return.max, plasma.only = T)
  LL.PARTNERSPREP = PARTNERSPREP.Likelihood(params, K, return.max = return.max, plasma.only = T)
  LL.VOICE + LL.FEMPREP + LL.PARTNERSPREP
}


FEMPREP.Likelihood = function(params, K, f = NULL, return.max = F, plasma.only = F, return.simulation = F){
  params$female = T
  
  FEMPREP.active.followup = 1024 + 1008 + 953 + 904 + 860 + 844 + 811 + 733 + 663 + 569 + 486 + 418 + 356 + 212
  FEMPREP.placebo.followup = 1032 + 1019 + 963 + 917 + 864 + 841 + 799 + 736 + 659 + 565 + 491 + 420 + 360 + 229
  FEMPREP.active.sero = 33
  FEMPREP.placebo.sero = 35

  
  FEMPREP.pos.control.samp = 35
  FEMPREP.neg.control.samp = 60
  FEMPREP.pos.sero = 7
  FEMPREP.neg.sero = 26
  
  
  positive.negative.plasma.Likelihood(params, K, female.pos.082.FEMPREP, female.neg.082.FEMPREP, f = f, plasma.only = plasma.only,
                                      placebo.sero = FEMPREP.placebo.sero, 
                                      placebo.followup = FEMPREP.placebo.followup,
                                      active.sero = FEMPREP.active.sero,
                                      active.followup = FEMPREP.active.followup,
                                      seroconverters.positive = FEMPREP.pos.sero,
                                      seroconverters.negative = FEMPREP.neg.sero,
                                      control.positive.sample = FEMPREP.pos.control.samp,
                                      control.negative.sample = FEMPREP.neg.control.samp,
                                      return.max = return.max,
                                      return.simulation = return.simulation)
}

VOICE.Likelihood = function(params, K, f = NULL, return.max = F, plasma.only = F, return.simulation = F){
  params$female = T
  
  VOICE.active.followup = 994 + 971 + 953 + 931 + 802 + 467 + 271 + 143 + 71 + 22 + 2
  VOICE.placebo.followup = 1008 + 982 +966 + 943 + 818 + 485 + 281 + 152 + 72 + 22 + 1
  VOICE.active.sero = 61
  VOICE.placebo.sero = 60
  
  VOICE.pos.control.samp = 77
  VOICE.neg.control.samp = 71
  VOICE.pos.sero = 24
  VOICE.neg.sero = 37
  
  positive.negative.plasma.Likelihood(params, K, female.pos.082.VOICE, female.neg.082.VOICE, f = f, plasma.only = plasma.only,
                                      placebo.sero = VOICE.placebo.sero, 
                                      placebo.followup = VOICE.placebo.followup,
                                      active.sero = VOICE.active.sero,
                                      active.followup = VOICE.active.followup,
                                      seroconverters.positive = VOICE.pos.sero,
                                      seroconverters.negative = VOICE.neg.sero,
                                      control.positive.sample = VOICE.pos.control.samp,
                                      control.negative.sample = VOICE.neg.control.samp,
                                      return.max = return.max,
                                      return.simulation = return.simulation)
}


PARTNERSPREP.Likelihood = function(params, K, f = NULL, return.max = F, plasma.only = F, return.simulation = F){
  params$female = T
  PARTNERSPREP.active.followup = round(4*9/.0095)
  PARTNERSPREP.placebo.followup = round(4*28/.0281)
  PARTNERSPREP.active.sero = 9
  PARTNERSPREP.placebo.sero = 28
  
  
  PARTNERSPREP.pos.control.samp = 135
  PARTNERSPREP.neg.control.samp = 40
  PARTNERSPREP.pos.sero = 1
  PARTNERSPREP.neg.sero = 7
  
  positive.negative.plasma.Likelihood(params, K, female.pos.082.PARTNERSPREP, female.neg.082.PARTNERSPREP, f = f, plasma.only = plasma.only, 
    placebo.sero = PARTNERSPREP.placebo.sero, 
    placebo.followup = PARTNERSPREP.placebo.followup,
    active.sero = PARTNERSPREP.active.sero,
    active.followup = PARTNERSPREP.active.followup,
    seroconverters.positive = PARTNERSPREP.pos.sero,
    seroconverters.negative = PARTNERSPREP.neg.sero,
    control.positive.sample = PARTNERSPREP.pos.control.samp,
    control.negative.sample = PARTNERSPREP.neg.control.samp,
    return.max = return.max,
    return.simulation = return.simulation)
}


positive.negative.plasma.Likelihood = function(params, K, positive.pk, negative.pk, f = NULL, plasma.only = F,
                                   placebo.sero = 0, 
                                   placebo.followup = 0,
                                   active.sero = 0,
                                   active.followup = 0,
                                   seroconverters.positive = 0,
                                   seroconverters.negative = 0,
                                   control.positive.sample = 0,
                                   control.negative.sample = 0,
                                   return.max = F,
                                   return.simulation = F){
  
  followups.positive = seq(seroconverters.positive, active.followup - seroconverters.negative)
  followups.negative = active.followup - followups.positive
  
  control.followup = active.followup - active.sero
  
  p.k.max = max(dbinom(seroconverters.positive, size = followups.positive, prob = seroconverters.positive/followups.positive, log = TRUE)
                + dbinom(seroconverters.negative, size = followups.negative, prob = seroconverters.negative/followups.negative, log = TRUE))
  
  if(return.max){
    return(p.k.max)
  }
  
  p=0
  for(i in seq(K)){
    
    #Randomly draw risk from beta distribution derived from placebo arm
    base.risk = rbeta(1, shape1 = placebo.sero + 1, shape2 = placebo.followup + 1 - placebo.sero)
    
    #Randomly draw positive plasma sample (for all members of uninfected active arm) from sub-sample
    control.positive.true = round(rbeta(1, shape1 = control.positive.sample + 1, shape2 =  control.negative.sample+1) * control.followup)
    
    positive.plasma = seroconverters.positive + control.positive.true
    negative.plasma = active.followup - positive.plasma
    
    positive.sample = 0
    negative.sample = 0
    
    if(plasma.only){
      positive.sample = params$positive.hazard * base.risk
      negative.sample = params$negative.hazard * base.risk
    }
    else{
      positive.sample = apply.hazard.ratio.func.resample(positive.pk, f, params, positive.plasma, base.risk)
      negative.sample = apply.hazard.ratio.func.resample(negative.pk, f, params, negative.plasma, base.risk)
    }
    
    if(return.simulation){
      return(c(positive = rbinom(1, size = positive.plasma, prob = positive.sample),
        negative = rbinom(1, size = negative.plasma, prob = negative.sample)))
    }
    
    p.k = dbinom(seroconverters.positive, size = positive.plasma, prob = positive.sample, log = TRUE) + dbinom(seroconverters.negative, size = negative.plasma, prob = negative.sample, log = TRUE)
    
    p = p + exp(p.k)/K
  }
  log(p)
}

# Helper Functions ####

sample.risk = function(risk, followup.interval){
  with(as.list(risk),{
    r = rnorm(1, mean = med, sd = se)
    1 - exp(-r*followup.interval)
  })
  
}


apply.hazard.ratio.func.resample = function(dbs, f, params, N, base.risk){
  dbs.N = sample(dbs, N, replace = T)
  HR = sapply(dbs.N, f, params)
  mean(1 - (1-base.risk)^HR, na.rm = T)
}

apply.hazard.ratio.func = function(dbs, f, params){
  mean(sapply(dbs, f, params))
}


