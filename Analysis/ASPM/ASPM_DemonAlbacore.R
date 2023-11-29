# ———————————————————————————————
# Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# ———————————————————————————————


########################################
## compile the model and mcmc options ##
########################################

ASPM <- cmdstanr::cmdstan_model('../../Models/ASPM/Main/main.stan',
                                pedantic=F, 
                                force_recompile=F)

seed=1234
nchains=nChains
nwarmup=nWarmupDemon
nsample=nSampleDemon
adaptDelta=0.99
maxTree=10

###########################################################################################################################
################################## Base case run (based on Punt et al., 1995) #############################################
###########################################################################################################################

###################
###### Prior ######
###################
PriorMu_transMh=c(log(M), bounded_logit(h, 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.5), 0.5) 

tau2Prior=c(1.708603, 0.008613854)
#tau2Prior=c(1/(0.5^2)+2,  0.0081*(1/(0.5^2)+2+1))
qPrior=c(1e-3, 1)

Sigma_MR0hRaw=diag(3)
#CorMR0=0.9
CorMR0=0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


R0Prior=c(exp(10), exp(20))


###############################
## Inputs and initial values ##
###############################

InputASPM=list("ntimes"=ntimes,
               "nages"=nages,
               "It"=It,
               "Ct"=Ct,
               "PriorMu_transMh"=PriorMu_transMh,
               "PriorSigma_transMh"=PriorSigma_transMh,
               "Sigma_MR0hRaw"=Sigma_MR0hRaw,
               "tau2Prior"=tau2Prior,
               "qPrior"=qPrior,
               "R0Prior"=R0Prior,
               "sigR"=sigR,
               "a_mat"=a_mat,
               "sig_mat"=sig_mat,
               "a_sel"=a_sel,
               "sig_sel"=sig_sel,
               "Linf"=Linf,
               "kappa"=kappa,
               "a0"=a0,
               "lwa"=lwa,
               "lwb"=lwb,
               "FemaleProp"=FemaleProp)


InitASPM=list(list("MR0hRaw"=c(0,0,0),
                   "logq"=log(1e-2),
                   "tau2"=0.01,
                   "logR0"=log(5e+6),
                   "logRt"=rep(log(1e+8), ntimes)))

##############
## MCMC run ##
##############

BaseModel=ASPM$sample(data=InputASPM, 
                      init = rep(InitASPM, nchains),
                      seed = seed,
                      chains = nchains, 
                      parallel_chains = nchains,  
                      iter_warmup = nwarmup, 
                      iter_sampling = nsample,
                      adapt_delta = adaptDelta,
                      max_treedepth = maxTree)


BaseModel=as_draws_matrix(BaseModel$draws())

#EstimASPM$summary(c("M","R0","h","q","tau2"))$median
