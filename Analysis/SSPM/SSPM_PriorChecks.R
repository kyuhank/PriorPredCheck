##########################################################################################################################################################
################################################### SBC (parameter values are drawn from a joint prior) ###################################################
##########################################################################################################################################################

QuantOfInterest=c("K","r", "q","sigma2" ,"tau2", "BStatus")

nPriorSamp=10000


#######################################
## true pars from the fitted results ##
#######################################

ParsTrueVals=c(268.2235, 0.292562, 0.236948, 0.002649985, 0.01145875)
names(ParsTrueVals)<-c("K","r","q","sigma2","tau2")
ParsTrueVals["sigma2"]=log(0.1^2+1)
ParsTrueVals["tau2"]=log(0.1^2+1)


##################################################################################################################################
########################################################## Prior (base) ##########################################################
##################################################################################################################################

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

BasePriors=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                "sigma2"=sigma2Prior,
                "tau2"=tau2Prior,
                "q"=qPrior) 

##################################################################################################################################
########################################################## Prior (biased) ########################################################
##################################################################################################################################

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

BiasedPriors=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
                  "sigma2"=sigma2Prior,
                  "tau2"=tau2Prior,
                  "q"=qPrior) 

###############################################################################################################################
########################################################## Prior (MVN) ########################################################
###############################################################################################################################

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

MVNPriors=list("logKr"=list(PriorMu_logKr, PriorSigma_logKr),
               "sigma2"=sigma2Prior,
               "tau2"=tau2Prior,
               "q"=qPrior) 


####################################################################################################################################
######################################################### Derive True quant ########################################################
####################################################################################################################################

ParsTrue=list("FixedQuant"=list("K"=ParsTrueVals["K"], 
                                "r"=ParsTrueVals["r"], 
                                "sigma2"=ParsTrueVals["sigma2"], 
                                "q"=ParsTrueVals["q"], 
                                "tau2"=ParsTrueVals["tau2"]))

TrueQuant=replicate(nPriorSamp, SspmSimData(ntimes=ntimes,
                                            Ct=Ct,
                                            ParsPrior = ParsTrue)[QuantOfInterest], simplify = F)

TrueQuant=do.call(rbind, lapply(TrueQuant, function(x) unlist(x[QuantOfInterest])))


###################################################################################################################################################
############################################################## Derive Effective Prior #############################################################
###################################################################################################################################################

#################
## Base Priors ##
#################

PriorSampBase=replicate(nPriorSamp, SspmSimData(ntimes=ntimes,
                                                Ct=Ct,
                                                ParsPrior = BasePriors)[QuantOfInterest], simplify = F)

PriorSampBase=do.call(rbind, lapply(PriorSampBase, function(x) unlist(x[QuantOfInterest])))


###################
## Biased Priors ##
###################

PriorSampBiased=replicate(nPriorSamp, SspmSimData(ntimes=ntimes,
                                                  Ct=Ct,
                                                  ParsPrior = BiasedPriors)[QuantOfInterest], simplify = F)

PriorSampBiased=do.call(rbind, lapply(PriorSampBiased, function(x) unlist(x[QuantOfInterest])))


################
## MVN Priors ##
################

PriorSampMVN=replicate(nPriorSamp, SspmSimData(ntimes=ntimes,
                                               Ct=Ct,
                                               ParsPrior = MVNPriors)[QuantOfInterest], simplify = F)

PriorSampMVN=do.call(rbind, lapply(PriorSampMVN, function(x) unlist(x[QuantOfInterest])))
