
##########################################################################################################################################################
################################################### SBC (parameter values are drawn from a joint prior) ###################################################
##########################################################################################################################################################

QuantOfInterest=c("M","R0","h", "q","tau2", "BStatus")

nPriorSamp=10000

#######################################
## true pars from the fitted results ##
#######################################

ParsTrueVals=c(0.2162035, 4263800, 0.849657, 0.233328, 0.01239915)
names(ParsTrueVals)<-c("M","R0","h","q","tau2")
ParsTrueVals["tau2"]=lnSigmaFromCV(0.1)^2

##################################################################################################################################
########################################################## Prior (base) ##########################################################
##################################################################################################################################

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


BasePriors=list("transMh"=list(PriorMu_transMh, PriorSigma_transMh),
                "tau2"=list(tau2Prior[1], tau2Prior[2]),
                "q"=list(qPrior[1], qPrior[2]),
                "R0"=list(R0Prior[1], R0Prior[2]),
                "Sigma_MR0hRaw"=Sigma_MR0hRaw)


###############################################################################################################################
###################################################  Prior (Fixed M)  #########################################################
###############################################################################################################################


PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.001), 0.5) 

tau2Prior=c(1/(0.5^2)+2,  0.00806*(1/(0.5^2)+2+1))
qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


FixedMPriors=list("transMh"=list(PriorMu_transMh, PriorSigma_transMh),
                  "tau2"=list(tau2Prior[1], tau2Prior[2]),
                  "q"=list(qPrior[1], qPrior[2]),
                  "R0"=list(R0Prior[1], R0Prior[2]),
                  "Sigma_MR0hRaw"=Sigma_MR0hRaw)



###############################################################################################################################
###################################################  Prior (Fixed R0)  ########################################################
###############################################################################################################################

PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(lnSigmaFromCV(0.5), 0.5) 

tau2Prior=c(1/(0.5^2)+2,  0.00806*(1/(0.5^2)+2+1))
qPrior=c(1e-3, 1)

R0CV=0.001

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0

FixedR0Priors=list("transMh"=list(PriorMu_transMh, PriorSigma_transMh),
                   "tau2"=list(tau2Prior[1], tau2Prior[2]),
                   "q"=list(qPrior[1], qPrior[2]),
                   "R0"=list(R0Prior[1], R0Prior[2]),
                   "Sigma_MR0hRaw"=Sigma_MR0hRaw)




###############################################################################################################################
###################################################  Prior (Biased)  #########################################################
###############################################################################################################################


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


BiasedPriors=list("transMh"=list(PriorMu_transMh, PriorSigma_transMh),
                  "tau2"=list(tau2Prior[1], tau2Prior[2]),
                  "q"=list(qPrior[1], qPrior[2]),
                  "R0"=list(R0Prior[1], R0Prior[2]),
                  "Sigma_MR0hRaw"=Sigma_MR0hRaw)


##########################################################################################################################################
###################################################  Prior (estimating R0 only)  #########################################################
##########################################################################################################################################


PriorMu_transMh=c(log(ParsTrueVals["M"]), bounded_logit(ParsTrueVals["h"], 0.2, 1))
PriorSigma_transMh=c(1e-3, 1e-3) 

tau2Prior=c(1/(0.5^2)+2,  0.00806*(1/(0.5^2)+2+1))
qPrior=c(1e-3, 1)

R0CV=0.2

width=(R0CV*sqrt(12))/2
R0Prior=c(exp(log(ParsTrueVals["R0"])*(1-width)), exp(log(ParsTrueVals["R0"])*(1+width)))

Sigma_MR0hRaw=diag(3)
CorMR0=0.0
Sigma_MR0hRaw[1,2]=CorMR0
Sigma_MR0hRaw[2,1]=CorMR0


R0onlyPriors=list("transMh"=list(PriorMu_transMh, PriorSigma_transMh),
                  "tau2"=list(tau2Prior[1], tau2Prior[2]),
                  "q"=list(qPrior[1], qPrior[2]),
                  "R0"=list(R0Prior[1], R0Prior[2]),
                  "Sigma_MR0hRaw"=Sigma_MR0hRaw)

####################################################################################################################################
######################################################### Derive True quant ########################################################
####################################################################################################################################

TrueQuant=replicate(nPriorSamp, AspmSimData(a_mat=a_mat, 
                                            nages=nages, 
                                            sig_mat=sig_mat, 
                                            kappa=kappa, 
                                            Linf=Linf, 
                                            a0=a0, 
                                            lwa=lwa,
                                            lwb=lwb, 
                                            a_sel=a_sel,
                                            sig_sel=sig_sel, 
                                            sigR=sigR, 
                                            Ct=Ct, 
                                            ObtainRef = F,
                                            ParsPrior = list("FixedQuant"=as.list(ParsTrueVals)))[QuantOfInterest], simplify = F)

TrueQuant=do.call(rbind, lapply(TrueQuant, function(x) unlist(x[QuantOfInterest])))

###################################################################################################################################################
############################################################## Derive Effective Prior #############################################################
###################################################################################################################################################

######################################
############# Base Priors ############
######################################

PriorSampBase=replicate(nPriorSamp, AspmSimData(a_mat=a_mat, 
                                                nages=nages, 
                                                sig_mat=sig_mat, 
                                                kappa=kappa, 
                                                Linf=Linf, 
                                                a0=a0, 
                                                lwa=lwa,
                                                lwb=lwb, 
                                                a_sel=a_sel,
                                                sig_sel=sig_sel, 
                                                sigR=sigR, 
                                                Ct=Ct, 
                                                ObtainRef = F,
                                                ParsPrior = BasePriors)[QuantOfInterest], simplify = F)

PriorSampBase=do.call(rbind, lapply(PriorSampBase, function(x) unlist(x[QuantOfInterest])))



########################################
############# FixedM Priors ############
########################################

PriorSampFixedM=replicate(nPriorSamp, AspmSimData(a_mat=a_mat, 
                                                  nages=nages, 
                                                  sig_mat=sig_mat, 
                                                  kappa=kappa, 
                                                  Linf=Linf, 
                                                  a0=a0, 
                                                  lwa=lwa,
                                                  lwb=lwb, 
                                                  a_sel=a_sel,
                                                  sig_sel=sig_sel, 
                                                  sigR=sigR, 
                                                  Ct=Ct, 
                                                  ObtainRef = F,
                                                  ParsPrior = FixedMPriors)[QuantOfInterest], simplify = F)

PriorSampFixedM=do.call(rbind, lapply(PriorSampFixedM, function(x) unlist(x[QuantOfInterest])))




########################################
############# FixedR0 Priors ############
########################################

PriorSampFixedR0=replicate(nPriorSamp, AspmSimData(a_mat=a_mat, 
                                                   nages=nages, 
                                                   sig_mat=sig_mat, 
                                                   kappa=kappa, 
                                                   Linf=Linf, 
                                                   a0=a0, 
                                                   lwa=lwa,
                                                   lwb=lwb, 
                                                   a_sel=a_sel,
                                                   sig_sel=sig_sel, 
                                                   sigR=sigR, 
                                                   Ct=Ct, 
                                                   ObtainRef = F,
                                                   ParsPrior = FixedR0Priors)[QuantOfInterest], simplify = F)

PriorSampFixedR0=do.call(rbind, lapply(PriorSampFixedR0, function(x) unlist(x[QuantOfInterest])))




########################################
############# Biased Priors ############
########################################

PriorSampBiased=replicate(nPriorSamp, AspmSimData(a_mat=a_mat, 
                                                  nages=nages, 
                                                  sig_mat=sig_mat, 
                                                  kappa=kappa, 
                                                  Linf=Linf, 
                                                  a0=a0, 
                                                  lwa=lwa,
                                                  lwb=lwb, 
                                                  a_sel=a_sel,
                                                  sig_sel=sig_sel, 
                                                  sigR=sigR, 
                                                  Ct=Ct, 
                                                  ObtainRef = F,
                                                  ParsPrior = BiasedPriors)[QuantOfInterest], simplify = F)

PriorSampBiased=do.call(rbind, lapply(PriorSampBiased, function(x) unlist(x[QuantOfInterest])))







########################################
############# R0only Priors ############
########################################

PriorSampR0only=replicate(nPriorSamp, AspmSimData(a_mat=a_mat, 
                                                  nages=nages, 
                                                  sig_mat=sig_mat, 
                                                  kappa=kappa, 
                                                  Linf=Linf, 
                                                  a0=a0, 
                                                  lwa=lwa,
                                                  lwb=lwb, 
                                                  a_sel=a_sel,
                                                  sig_sel=sig_sel, 
                                                  sigR=sigR, 
                                                  Ct=Ct, 
                                                  ObtainRef = F,
                                                  ParsPrior = R0onlyPriors)[QuantOfInterest], simplify = F)

PriorSampR0only=do.call(rbind, lapply(PriorSampR0only, function(x) unlist(x[QuantOfInterest])))

