

# ———————————————————————————————
# Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# ———————————————————————————————


########################################
## compile the model and mcmc options ##
########################################

SSPM <- cmdstanr::cmdstan_model('../../src/SSPM/SSPM.stan',
                                pedantic=F, 
                                force_recompile=F)

###################################################################################################################################
###################################################  Parametric bootstrap test  ###################################################
###################################################################################################################################

seed=12345
nchains=nChains
nwarmup=nWarmup
nsample=nSample
adaptDelta=0.99
maxTree=10

nsim=nSimBoot

QuantOfInterest=c("K","r","q","sigma2", "tau2", "BStatus", "Bt", "Ht")


Inits=list(list("KrRaw"=c(1, 1),
                "logq"=log(0.01),
                "sigma2"=0.1,
                "tau2"=0.1,
                "Pt"=rep(1, ntimes+1),
                "Ht"=rep(0.05, ntimes)))


#######################################
## true pars from the fitted results ##
#######################################

#ParsTrueVals=BaseModel$summary(c("K","r","q","sigma2","tau2"))$median
ParsTrueVals=c(268.2235, 0.292562, 0.236948, 0.002649985, 0.01145875)
names(ParsTrueVals)<-c("K","r","q","sigma2","tau2")
ParsTrueVals["sigma2"]=log(0.1^2+1)
ParsTrueVals["tau2"]=log(0.1^2+1)

###################################################################################################################
###################################################  Base case  ###################################################
###################################################################################################################

##########################
###### Prior (base) ######
##########################

CVr=0.5
CVk=0.5

## prior for K and r

PriorMu_logKr=c(log(ParsTrueVals["K"]), log(ParsTrueVals["r"]))
SDs=diag(c(sqrt(log(CVk^2+1)), sqrt(log(CVr^2+1))))


CoRs=diag(2)
CorKr=0
CoRs[1,2]=CorKr
CoRs[2,1]=CorKr

PriorSigma_logKr=SDs%*%CoRs%*%SDs


## prior for variance parameters

CVtau=0.5
CVsig=0.5

tau2Prior=c(1/(CVtau^2)+2,  0.00806*(1/(CVtau^2)+2+1))
sigma2Prior=c(1/(CVsig^2)+2,  0.00806*(1/(CVsig^2)+2+1))

qPrior=c(1e-3, 1)

## inputs for the base model

InputBaseBoot=list("ntimes"=ntimes,
                   "It"=It,
                   "Ct"=Ct,
                   "PriorMu_logKr"=PriorMu_logKr,
                   "PriorSigma_logKr"=PriorSigma_logKr,
                   "sigma2Prior"=sigma2Prior,
                   "tau2Prior"=tau2Prior,
                   "qPrior"=qPrior)

BasePriorsBoot=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                    "sigma2"=sigma2Prior,
                    "tau2"=tau2Prior,
                    "q"=qPrior) 

###################
## Bootstrap run ##
###################

SSPM_Boot_base=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                   ModelType="SSPM",
                   StanObj=SSPM,
                   QuantOfInterest=QuantOfInterest,
                   InputData=InputBaseBoot, 
                   InitPars=Inits,
                   seed=seed,
                   nchains=nchains,
                   nwarmup=nwarmup,
                   nsample=nsample,
                   adaptDelta=adaptDelta,
                   maxTree=maxTree,
                   nsim=nsim)


###################################################################################################################
###################################################  Biased case  #################################################
###################################################################################################################

############################
###### Prior (biased) ######
############################

CVr=0.5
CVk=0.5

## prior for K and r

PriorMu_logKr=c(log(ParsTrueVals["K"]-150), log(ParsTrueVals["r"]-0.1))
SDs=diag(c(sqrt(log(CVk^2+1)), sqrt(log(CVr^2+1))))


CoRs=diag(2)
CorKr=0
CoRs[1,2]=CorKr
CoRs[2,1]=CorKr

PriorSigma_logKr=SDs%*%CoRs%*%SDs


## prior for variance parameters

CVtau=0.5
CVsig=0.5

tau2Prior=c(1/(CVtau^2)+2,  0.00806*(1/(CVtau^2)+2+1))
sigma2Prior=c(1/(CVsig^2)+2,  0.00806*(1/(CVsig^2)+2+1))

qPrior=c(1e-3, 1)

## inputs for the base model

InputBiasedBoot=list("ntimes"=ntimes,
                     "It"=It,
                     "Ct"=Ct,
                     "PriorMu_logKr"=PriorMu_logKr,
                     "PriorSigma_logKr"=PriorSigma_logKr,
                     "sigma2Prior"=sigma2Prior,
                     "tau2Prior"=tau2Prior,
                     "qPrior"=qPrior)

BiasedPriorsBoot=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                      "sigma2"=sigma2Prior,
                      "tau2"=tau2Prior,
                      "q"=qPrior) 

###################
## Bootstrap run ##
###################

SSPM_Boot_biased=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                     ModelType="SSPM",
                     StanObj=SSPM,
                     QuantOfInterest=QuantOfInterest,
                     InputData=InputBiasedBoot, 
                     InitPars=Inits,
                     seed=seed,
                     nchains=nchains,
                     nwarmup=nwarmup,
                     nsample=nsample,
                     adaptDelta=adaptDelta,
                     maxTree=maxTree,
                     nsim=nsim)

###################################################################################################################
###################################################  MVN case  #################################################
###################################################################################################################

############################
###### Prior (biased) ######
############################

CVr=0.5
CVk=0.3

## prior for K and r

PriorMu_logKr=c(log(ParsTrueVals["K"]), log(ParsTrueVals["r"]))
SDs=diag(c(sqrt(log(CVk^2+1)), sqrt(log(CVr^2+1))))


CoRs=diag(2)
CorKr=-0.9
CoRs[1,2]=CorKr
CoRs[2,1]=CorKr

PriorSigma_logKr=SDs%*%CoRs%*%SDs


## prior for variance parameters

CVtau=0.5
CVsig=0.5

tau2Prior=c(1/(CVtau^2)+2,  0.00806*(1/(CVtau^2)+2+1))
sigma2Prior=c(1/(CVsig^2)+2,  0.00806*(1/(CVsig^2)+2+1))

qPrior=c(1e-3, 1)

## inputs for the base model

InputMVNBoot=list("ntimes"=ntimes,
                   "It"=It,
                   "Ct"=Ct,
                   "PriorMu_logKr"=PriorMu_logKr,
                   "PriorSigma_logKr"=PriorSigma_logKr,
                   "sigma2Prior"=sigma2Prior,
                   "tau2Prior"=tau2Prior,
                   "qPrior"=qPrior)

MVNPriorsBoot=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                   "sigma2"=sigma2Prior,
                   "tau2"=tau2Prior,
                   "q"=qPrior) 

###################
## Bootstrap run ##
###################

SSPM_Boot_MVN=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                  ModelType="SSPM",
                  StanObj=SSPM,
                  QuantOfInterest=QuantOfInterest,
                  InputData=InputMVNBoot, 
                  InitPars=Inits,
                  seed=seed,
                  nchains=nchains,
                  nwarmup=nwarmup,
                  nsample=nsample,
                  adaptDelta=adaptDelta,
                  maxTree=maxTree,
                  nsim=nsim)


############################
########## save ############
############################

save.image(file="SSPM_analysis.RData")

