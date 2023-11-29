
# ———————————————————————————————
# Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# ———————————————————————————————


########################################
## compile the model and mcmc options ##
########################################

SSPM <- cmdstanr::cmdstan_model('../../Models/SSPM/SSPM.stan',
                                pedantic=F, 
                                force_recompile=F)

seed=12345
nchains=nChains
nwarmup=nWarmupDemon
nsample=nSampleDemon
adaptDelta=0.999
maxTree=10

################################################################################################################################
################################## Base case run (based on Millar and Meyer, 2000) #############################################
################################################################################################################################

Inits=list(list("KrRaw"=c(1, 1),
                "logq"=log(0.01),
                "sigma2"=0.01,
                "tau2"=0.01,
                "Pt"=rep(1, ntimes+1),
                "Ht"=rep(0.05, ntimes)))


################################################
###### Prior (from Millar and Meyer, 2000) #####
################################################

PriorMu_logKr=c(5.042905, -1.38)
PriorSigma_logKr=diag(c(1/sqrt(3.7603664), 1/sqrt(3.845))^2)

sigma2Prior=c(3.785518, 0.010223)
tau2Prior=c(1.708603, 0.008613854)
qPrior=c(1e-3, 1)


BasePriors=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                "sigma2"=sigma2Prior,
                "tau2"=tau2Prior,
                "q"=qPrior)

###############################
## Inputs and initial values ##
###############################

InputBase=list("ntimes"=ntimes,
               "It"=It,
               "Ct"=Ct,
               "PriorMu_logKr"=PriorMu_logKr,
               "PriorSigma_logKr"=PriorSigma_logKr,
               "sigma2Prior"=sigma2Prior,
               "tau2Prior"=tau2Prior,
               "qPrior"=qPrior)

##############
## MCMC run ##
##############

BaseModel=SSPM$sample(data=InputBase, 
                      init = rep(Inits, nchains),
                      seed = seed,  
                      #output_dir = CSVdir, 
                      #output_basename ="BaseModel",
                      chains = nchains, 
                      parallel_chains = nchains,  
                      iter_warmup = nwarmup, 
                      iter_sampling = nsample,
                      adapt_delta = adaptDelta,
                      max_treedepth = maxTree)


BasePost=as_draws_matrix(BaseModel$draws())

##################################################################################################################################################
###################################################### Find Effective Prior ######################################################################
##################################################################################################################################################

nPriorSamp=10000

BaseInputPriors=cbind("K"=rlnorm(nPriorSamp, BasePriors$logKr[[1]][1], BasePriors$logKr[[2]][1,1]),
                      "r"=rlnorm(nPriorSamp, BasePriors$logKr[[1]][2], BasePriors$logKr[[2]][2,2]),
                      "sigma2"=1/rgamma(nPriorSamp, sigma2Prior[1], sigma2Prior[2]),
                      "q"=exp(runif(nPriorSamp, log(BasePriors$q[1]), log(BasePriors$q[1]))),
                      "tau2"=1/rgamma(nPriorSamp, sigma2Prior[1], sigma2Prior[2]))

BaseEffPriors=do.call(rbind, replicate(nPriorSamp, unlist(SspmSimData(BasePriors,
                                                                      ntimes,
                                                                      Ct)[c("K","r","sigma2") ]), simplify = F))


#############################################################
######## approximate the effective prior using MVN ##########
#############################################################

PriorMu_logKr_mvn=apply(log(BaseEffPriors[,c("K","r")] ), 2, mean)
PriorSigma_logKr_mvn=cov(log(BaseEffPriors[,c("K","r")]))

InputMVN=list("ntimes"=ntimes,
              "It"=It,
              "Ct"=Ct,
              "PriorMu_logKr"=PriorMu_logKr_mvn,
              "PriorSigma_logKr"=PriorSigma_logKr_mvn,
              "sigma2Prior"=sigma2Prior,
              "tau2Prior"=tau2Prior,
              "qPrior"=qPrior)

##############
## MCMC run ##
##############

MvnModel=SSPM$sample(data=InputMVN, 
                     init = rep(Inits, nchains),
                     seed = seed, 
                     chains = nchains, 
                     #output_dir = CSVdir,
                     #output_basename ="MvnModel",
                     parallel_chains = nchains,  
                     iter_warmup = nwarmup, 
                     iter_sampling = nsample,
                     adapt_delta = adaptDelta,
                     max_treedepth = maxTree)


MvnPost=as_draws_matrix(MvnModel$draws())


###################################################################################################################################################
###################################################### Reparameterised version ####################################################################
###################################################################################################################################################

#################################################################################################################################
############################################### Model where beta(1, 1) prior on Ht ##############################################
#################################################################################################################################

########################################
## compile the model and mcmc options ##
########################################

SSPM_explicit_flat <- cmdstanr::cmdstan_model('../../Models/SSPM/SSPM_explicit_Flat.stan',
                                              pedantic=F, 
                                              force_recompile=F)


Inits_explicit_flat=list(list("KrRaw"=c(1, 1),
                              "logq"=log(0.01),
                              "sigma2"=0.01,
                              "tau2"=0.01,
                              "Pt"=rep(1, ntimes+1),
                              "Ht"=rep(0.05, ntimes)))

################################################
###### Prior (from Millar and Meyer, 2000) #####
################################################

PriorMu_logKr=c(5.042905, -1.38)
PriorSigma_logKr=diag(c(1/sqrt(3.7603664), 1/sqrt(3.845))^2)

sigma2Prior=c(3.785518, 0.010223)
tau2Prior=c(1.708603, 0.008613854)
qPrior=c(1e-3, 1)
HPrior=rbind(rep(1,ntimes), rep(1, ntimes))


Priors_explicit_flat=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                          "sigma2"=sigma2Prior,
                          "tau2"=tau2Prior,
                          "q"=qPrior,
                          "HPrior"=HPrior)

###############################
## Inputs and initial values ##
###############################

Input_explicit_flat=list("ntimes"=ntimes,
                         "It"=It,
                         "Ct"=Ct,
                         "PriorMu_logKr"=PriorMu_logKr,
                         "PriorSigma_logKr"=PriorSigma_logKr,
                         "sigma2Prior"=sigma2Prior,
                         "tau2Prior"=tau2Prior,
                         "qPrior"=qPrior,
                         "HPrior"=HPrior)

##############
## MCMC run ##
##############

ExplicitModel_flat=SSPM_explicit_flat$sample(data=Input_explicit_flat, 
                                             init = rep(Inits_explicit_flat, nchains),
                                             seed = seed,  
                                             #output_dir = CSVdir,
                                             #output_basename ="ExplicitModel_flat",
                                             chains = nchains, 
                                             parallel_chains = nchains,  
                                             iter_warmup = nwarmup, 
                                             iter_sampling = nsample,
                                             adapt_delta = adaptDelta,
                                             max_treedepth = maxTree+2)


capture.output({SSPM_explicit_flat_gen=SSPM_explicit_flat$generate_quantities(ExplicitModel_flat,
                                                                              data=Input_explicit_flat,
                                                                              parallel_chains = nchains,
                                                                              #output_dir = CSVdir, 
                                                                              #output_basename = "SSPM_explicit_flat_gen"
                                                                              )}, file=nullfile())

SSPM_explicit_flat_gen=as_draws_matrix(SSPM_explicit_flat_gen$draws())

ExplicitModel_flatPost=as_draws_matrix(ExplicitModel_flat$draws())


#################################################################################################################################
############################################### Model where beta(1, 7) prior on Ht ##############################################
#################################################################################################################################

########################################
## compile the model and mcmc options ##
########################################

SSPM_explicit_inform <- cmdstanr::cmdstan_model('../../Models/SSPM/SSPM_explicit_Flat.stan',
                                                pedantic=F, 
                                                force_recompile=F)


Inits_explicit_inform=list(list("KrRaw"=c(1, 1),
                              "logq"=log(0.01),
                              "sigma2"=0.01,
                              "tau2"=0.01,
                              "Pt"=rep(1, ntimes+1),
                              "Ht"=rep(0.05, ntimes)))

################################################
###### Prior (from Millar and Meyer, 2000) #####
################################################

PriorMu_logKr=c(5.042905, -1.38)
PriorSigma_logKr=diag(c(1/sqrt(3.7603664), 1/sqrt(3.845))^2)

sigma2Prior=c(3.785518, 0.010223)
tau2Prior=c(1.708603, 0.008613854)
qPrior=c(1e-3, 1)
HPrior=rbind(rep(1,ntimes), rep(7, ntimes))


Priors_explicit_inform=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                          "sigma2"=sigma2Prior,
                          "tau2"=tau2Prior,
                          "q"=qPrior,
                          "HPrior"=HPrior)

###############################
## Inputs and initial values ##
###############################

Input_explicit_inform=list("ntimes"=ntimes,
                         "It"=It,
                         "Ct"=Ct,
                         "PriorMu_logKr"=PriorMu_logKr,
                         "PriorSigma_logKr"=PriorSigma_logKr,
                         "sigma2Prior"=sigma2Prior,
                         "tau2Prior"=tau2Prior,
                         "qPrior"=qPrior,
                         "HPrior"=HPrior)

##############
## MCMC run ##
##############

ExplicitModel_inform=SSPM_explicit_inform$sample(data=Input_explicit_inform, 
                                             init = rep(Inits_explicit_inform, nchains),
                                             seed = seed,  
                                             #output_dir = CSVdir,
                                             #output_basename ="ExplicitModel_inform",
                                             chains = nchains, 
                                             parallel_chains = nchains,  
                                             iter_warmup = nwarmup, 
                                             iter_sampling = nsample,
                                             adapt_delta = adaptDelta,
                                             max_treedepth = maxTree+2)


capture.output({SSPM_explicit_inform_gen=SSPM_explicit_inform$generate_quantities(ExplicitModel_inform, 
                                                                                  data=Input_explicit_inform, 
                                                                                  parallel_chains = nchains, 
                                                                                  #output_dir = CSVdir, 
                                                                                  #output_basename = "SSPM_explicit_inform_gen"
                                                                                  )}, file=nullfile())

SSPM_explicit_inform_gen=as_draws_matrix(SSPM_explicit_inform_gen$draws())


ExplicitModel_informPost=as_draws_matrix(ExplicitModel_inform$draws())

