
########################################
## compile the model and mcmc options ##
########################################

SSPM <- cmdstanr::cmdstan_model('../../Models/SSPM/SSPM.stan',
                                pedantic=F, 
                                force_recompile=F)

##########################################################################################################################################################
################################################### SBC (parameter values are drawn from a joint prior) ###################################################
##########################################################################################################################################################
seed=12345
nchains=nChains
nwarmup=nWarmup
nsample=nSample
adaptDelta=0.999
maxTree=10

nsim=nSimSBC
QuantOfInterest=c("K","r","q","sigma2", "tau2", "BStatus", "Bt", "Ht")

###############################
## Inputs and initial values ##
###############################

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


Inits=list(list("KrRaw"=c(1, 1),
                "logq"=log(0.01),
                "sigma2"=0.01,
                "tau2"=0.01,
                "Pt"=rep(1, ntimes+1),
                "Ht"=rep(0.05, ntimes)))


###############################
############ SBC run ##########
###############################

SSPM_SBC=SBC(ParsPrior=BasePriors,
             ModelType="SSPM",
             StanObj=SSPM,
             QuantOfInterest=QuantOfInterest,
             InputData=InputBase, 
             InitPars=Inits,
             seed=seed,
             nchains=nchains,
             nwarmup=nwarmup,
             nsample=nsample,
             adaptDelta=adaptDelta,
             maxTree=maxTree,
             nPost = 10,
             nsim=nsim)

