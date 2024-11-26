
# ———————————————————————————————
# Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# ———————————————————————————————

########################################
## compile the model and mcmc options ##
########################################

ASPM <- cmdstanr::cmdstan_model('../../src/ASPM/Main/main.stan',
                                pedantic=F, 
                                force_recompile=F)

###################################################################################################################################
###################################################  Parametric bootstrap test  ###################################################
###################################################################################################################################

seed=123
nchains=nChains
nwarmup=nWarmup
nsample=nSample
adaptDelta=0.99
maxTree=10

nsim=nSimBoot
QuantOfInterest=c("M","R0","h", "q","tau2", "BStatus", "SSBStatus", "availBt", "Ht")

#######################################
## true pars from the fitted results ##
#######################################

ParsTrueVals=c(0.2162035, 4263800, 0.849657, 0.233328, 0.01239915)
names(ParsTrueVals)<-c("M","R0","h","q","tau2")

Inits=list(list("MR0hRaw"=c(0,0,0),
                "logq"=log(1e-2),
                "tau2"=0.01,
                "logR0"=log(5e+6),
                "logRt"=rep(log(1e+8), ntimes)))


##################################################################################################################################
###################################################  Matched CV=0.1 (high info) ##################################################
##################################################################################################################################


cat('\n\n Run Matched CV=0.1 (high info) \n\n')  

ParsTrueVals["tau2"]= lnSigmaFromCV(0.1)^2

####################################
###### Prior (Matched CV=0.1) ######
####################################

PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.5), 0.5) 

tau2Prior=c(1/(0.5^2)+2,  0.00806*(1/(0.5^2)+2+1))

qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


###############################
## Inputs and initial values ##
###############################

InputMatchedCV0.1Boot=list("ntimes"=ntimes,
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

###################
## Bootstrap run ##
###################

ASPM_Boot_MatchedCV0.1=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                           ModelType="ASPM",
                           StanObj=ASPM,
                           QuantOfInterest=QuantOfInterest,
                           InputData=InputMatchedCV0.1Boot, 
                           InitPars=Inits,
                           seed=seed,
                           nchains=nchains,
                           nwarmup=nwarmup,
                           nsample=nsample,
                           adaptDelta=adaptDelta,
                           maxTree=maxTree,
                           nsim=nsim)

##################################################################################################################################
###################################################  Matched CV=0.3 (low info) ###################################################
##################################################################################################################################

cat('\n\n Run Matched CV=0.3 (low info) \n\n')  

ParsTrueVals["tau2"]= lnSigmaFromCV(0.3)^2

####################################
###### Prior (Matched CV=0.1) ######
####################################

PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.5), 0.5) 

tau2Prior=c(1/(0.5^2)+2,  0.0698*(1/(0.5^2)+2+1))

qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


###############################
## Inputs and initial values ##
###############################

InputMatchedCV0.3Boot=list("ntimes"=ntimes,
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

###################
## Bootstrap run ##
###################

ASPM_Boot_MatchedCV0.3=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                           ModelType="ASPM",
                           StanObj=ASPM,
                           QuantOfInterest=QuantOfInterest,
                           InputData=InputMatchedCV0.3Boot, 
                           InitPars=Inits,
                           seed=seed,
                           nchains=nchains,
                           nwarmup=nwarmup,
                           nsample=nsample,
                           adaptDelta=adaptDelta,
                           maxTree=maxTree,
                           nsim=nsim)

######################################################################################################################################
###################################################  Biased CV=0.1 (high info case)  #################################################
######################################################################################################################################

cat('\n\n Run biased CV=0.1 (high info) \n\n')  

ParsTrueVals["tau2"]= lnSigmaFromCV(0.1)^2

###################################
###### Prior (Biased CV=0.1) ######
###################################

PriorMu_transMh=c(log(ParsTrueVals["M"]+0.15), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.5), 0.5) 

tau2Prior=c(1/(0.5^2)+2,  0.00806*(1/(0.5^2)+2+1))
qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)-5 ))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


###############################
## Inputs and initial values ##
###############################

InputBiasedCV0.1Boot=list("ntimes"=ntimes,
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

###################
## Bootstrap run ##
###################

ASPM_Boot_BiasedCV0.1=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                          ModelType="ASPM",
                          StanObj=ASPM,
                          QuantOfInterest=QuantOfInterest,
                          InputData=InputBiasedCV0.1Boot, 
                          InitPars=Inits,
                          seed=seed,
                          nchains=nchains,
                          nwarmup=nwarmup,
                          nsample=nsample,
                          adaptDelta=adaptDelta,
                          maxTree=maxTree,
                          nsim=nsim)

######################################################################################################################################
###################################################  Biased CV=0.3 (low info case)  #################################################
######################################################################################################################################

cat('\n\n Run biased CV=0.3 (low info) \n\n')  

ParsTrueVals["tau2"]= lnSigmaFromCV(0.3)^2

###################################
###### Prior (Biased CV=0.3) ######
###################################

PriorMu_transMh=c(log(ParsTrueVals["M"]+0.15), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.5), 0.5) 

tau2Prior=c(1/(0.5^2)+2,  0.0698*(1/(0.5^2)+2+1))
qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)-5 ))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


###############################
## Inputs and initial values ##
###############################

InputBiasedCV0.3Boot=list("ntimes"=ntimes,
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

###################
## Bootstrap run ##
###################

ASPM_Boot_BiasedCV0.3=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                          ModelType="ASPM",
                          StanObj=ASPM,
                          QuantOfInterest=QuantOfInterest,
                          InputData=InputBiasedCV0.3Boot, 
                          InitPars=Inits,
                          seed=seed,
                          nchains=nchains,
                          nwarmup=nwarmup,
                          nsample=nsample,
                          adaptDelta=adaptDelta,
                          maxTree=maxTree,
                          nsim=nsim)


##################################################################################################################################
###################################################  OnlyR0 CV=0.1 (high info) ##################################################
##################################################################################################################################

cat('\n\n Run OnlyR0 CV=0.1 (high info) \n\n')  

ParsTrueVals["tau2"]= lnSigmaFromCV(0.1)^2

####################################
###### Prior (OnlyR0 CV=0.1) #######
####################################

PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(1e-3), 1e-3) 

tau2Prior=c(1/(0.5^2)+2,  0.00806*(1/(0.5^2)+2+1))

qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


###############################
## Inputs and initial values ##
###############################

InputOnlyR0CV0.1Boot=list("ntimes"=ntimes,
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

###################
## Bootstrap run ##
###################

ASPM_Boot_OnlyR0CV0.1=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                          ModelType="ASPM",
                          StanObj=ASPM,
                          QuantOfInterest=QuantOfInterest,
                          InputData=InputOnlyR0CV0.1Boot, 
                          InitPars=Inits,
                          seed=seed,
                          nchains=nchains,
                          nwarmup=nwarmup,
                          nsample=nsample,
                          adaptDelta=adaptDelta,
                          maxTree=maxTree,
                          nsim=nsim)







##################################################################################################################################
###################################################  OnlyR0 CV=0.3 (low info) ##################################################
##################################################################################################################################

cat('\n\n Run OnlyR0 CV=0.3 (low info) \n\n')  

ParsTrueVals["tau2"]= lnSigmaFromCV(0.3)^2

####################################
###### Prior (OnlyR0 CV=0.3) #######
####################################

PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(1e-3), 1e-3) 

tau2Prior=c(1/(0.5^2)+2,  0.0698*(1/(0.5^2)+2+1))

qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


###############################
## Inputs and initial values ##
###############################

InputOnlyR0CV0.3Boot=list("ntimes"=ntimes,
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

###################
## Bootstrap run ##
###################

ASPM_Boot_OnlyR0CV0.3=SBC(ParsPrior=list("FixedQuant"=as.list(ParsTrueVals)),
                          ModelType="ASPM",
                          StanObj=ASPM,
                          QuantOfInterest=QuantOfInterest,
                          InputData=InputOnlyR0CV0.3Boot, 
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

save.image(file="ASPM_analysis.RData")