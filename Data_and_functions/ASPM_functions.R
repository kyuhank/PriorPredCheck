

########################
## selectivity-at-age ##
########################

selectivity_at_age=function (a_sel, nages, sig_sel) {
  
  a=1:nages
  
  Sel=1/(1+exp( (-(a-a_sel)/sig_sel) ))
  
  return("sel_a"=Sel)
  
}

#####################
## maturity-at-age ##
#####################

maturity_at_age=function (a_mat, nages, sig_mat) {
  
  #Mat=1/(1+exp( (-(a-a_mat)/sig_mat) ))
  
  Mat=c()
  
  for(a in 1:nages) {
  if(a<a_mat) {
    Mat[a]=0
  } else if (a>=a_mat) { 
    Mat[a]=1
  }  
  }
    
  return("matu_a"=Mat)
  
}

###################
## weight-at-age ##
###################

weight_at_age=function (lwa, lwb, length_a) {
  
  wa=lwa*length_a^(lwb)
  
  return("weight_a"=wa)
  
}

###################
## length-at-age ##
###################

length_at_age=function (kappa, Linf, a0, nages) {
  
  a=1:nages
  
  La=Linf*(1-exp(-kappa*(a-a0)))
  
  return("length_a"=La)
  
}


########################
## stock recruitment  ##
########################
StockRecruitment= function (h,
                            R0,
                            SSB0,
                            SSBt) {
  
  R=(4.*R0*h*SSBt)/((1.-h)*SSB0+(5.*h-1.)*SSBt)
  
  return (R)
  
}


#############
## Eq surv ##
#############

eqSurvRate=function(M, 
                    nages) {
  
  SurAtAge=c()
  
  for (a in 1:nages) {
    if(a==1) {
      SurAtAge[a]=1
    } else if (a>1 && a<nages) {
      SurAtAge[a]=SurAtAge[a-1]*exp(-M)
    } else if (a==nages) {
      SurAtAge[a]=(SurAtAge[a-1]*exp(-M))/(1 -exp(-M))
    } 
    
  }
  
  return (SurAtAge)
  
}



############################
#### Ref calculation #######
############################


NumericMSY=function(nages,
                    R0,
                    h,
                    eqSurv,
                    Weight_a,
                    Maturity_a,
                    FemaleProp,
                    Selex_a,
                    M) {
  
  
  ntimes=200
  eqH=c()
  
  Nat=matrix(0., nrow=nages, ncol=ntimes+1)
  Cat=matrix(0., nrow=nages, ncol=ntimes)
  Bt=c()
  availBt=c()
  SSBt=c()
  BStatus=c()
  Ht=c()
  eqNa=R0*eqSurv
  Ct=matrix(0, nrow=ntimes, ncol=200)
  
  eqH[1]=0.00
  for(i in 2:200) {
    eqH[i]=eqH[i-1]+0.005
  }
  
  ##@@@@@@@@@@@@@@@@@@@
  ##@@ Eq quantities @@
  ##@@@@@@@@@@@@@@@@@@@
  
  SSB0=sum(FemaleProp*eqNa*Maturity_a*Weight_a)
  B0=sum(eqNa*Weight_a)
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ##@@ initial N (i.e., t=1) @@
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  Nat[,1]=eqNa;
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@
  ##@@ main loop (from t=2)@@
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@   
  
  for (i in 1:200) {
    for(t in 2:(ntimes+1) ) {
      SSBt[t-1]=sum(FemaleProp*Nat[,t-1]*Maturity_a*Weight_a)
      for (a in 1:nages) {
        if (a==1) {
          Nat[a,t]=StockRecruitment(h, R0, SSB0, SSBt[t-1])
        } else if (a>1 && a<nages) {
          Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*eqH[i])
        } else if (a==nages) {
          Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*eqH[i])+Nat[a, t-1]*exp(-M)*(1.0-Selex_a[a]*eqH[i])
        }
        Cat[a,t-1]=Nat[a,t-1]*Selex_a[a]*eqH[i]*Weight_a[a]
        if(is.null(Cat[a,t-1])) Cat[a,t-1]=0
      }
      Ct[t-1,i]=sum(Cat[,t-1])    
    }
  }
  
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ##@@@@ returning quantities @@@@
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  MSY=max(Ct[ntimes,])
  Hmsy=eqH[which.max(Ct[ntimes,])]
  Bmsy=MSY/Hmsy
  
  return (list("MSY"=MSY,
               "Hmsy"=Hmsy,
               "Bmsy"=Bmsy))
}


###################################
######### ASPM dynamics ###########
###################################


ASPM_dynamics=function(Maturity_a,
                       Weight_a,
                       Selex_a,
                       FemaleProp=0.5,
                       sigR,
                       nages,
                       Ct,
                       ParsPrior,
                       ObtainRef=T) {
  
  ntimes=length(Ct)
  Nat=matrix(0, nrow=nages, ncol=ntimes+1)
  Bt=c()
  availBt=c()
  SSBt=c()
  BStatus=c()
  Ht=c()

  
  repeat{
    
    Rdevs=rnorm(ntimes, 0, sigR)  
    
    if(is.null(ParsPrior[["FixedQuant"]] )) {  
      mvnMR0hRaw=MASS::mvrnorm(1, rep(0, 3), ParsPrior[["Sigma_MR0hRaw"]])
      
      M=exp( mvnMR0hRaw[1]*ParsPrior[["transMh"]][[2]][1]+ParsPrior[["transMh"]][[1]][1])
      R0=exp( pnorm(mvnMR0hRaw[2])*(log(ParsPrior[["R0"]][[2]])-log(ParsPrior[["R0"]][[1]]))+log(ParsPrior[["R0"]][[1]]))
      h=bounded_inv_logit( mvnMR0hRaw[3]*ParsPrior[["transMh"]][[2]][2]+ParsPrior[["transMh"]][[1]][2], Lower=0.2, Upper=1)
      
      q=exp(runif(1, log(ParsPrior[["q"]][[1]]), log(ParsPrior[["q"]][[2]])) )
      tau2=1/rgamma(1, ParsPrior[["tau2"]][[1]], ParsPrior[["tau2"]][[2]])
      
    } else {
      
      M=ParsPrior[["FixedQuant"]][["M"]]  
      R0=ParsPrior[["FixedQuant"]][["R0"]]  
      h=ParsPrior[["FixedQuant"]][["h"]]
      q=ParsPrior[["FixedQuant"]][["q"]]  
      tau2=ParsPrior[["FixedQuant"]][["tau2"]]  
      
    }
    
    eqSurv=eqSurvRate(M, nages)
    
    eqNa=R0*eqSurv
    
    
    availB0=R0*sum(eqSurv*Weight_a*Selex_a)
    
    ##@@@@@@@@@@@@@@@@@@@
    ##@@ Eq quantities @@
    ##@@@@@@@@@@@@@@@@@@@
    
    SSB0=sum(FemaleProp*eqNa*Maturity_a*Weight_a)
    B0=sum(eqNa*Weight_a)
    
    Nat[,1]=eqNa
    
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@
    ##@@ main loop (from t=2)@@
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@   
    
    for(t in 2:(ntimes+1) ) {
      
      SSBt[t-1]=sum(FemaleProp*Nat[,t-1]*Maturity_a*Weight_a)
      availBt[t-1]=sum(Nat[,t-1]*Weight_a*Selex_a)
      Bt[t-1]=sum(Nat[,t-1]*Weight_a)
      
      Ht[t-1]=Ct[t-1]/availBt[t-1]
      
      for (a in 1:nages) {
        if (a==1) {
          Nat[a,t]=StockRecruitment(h, R0, SSB0, SSBt[t-1])*exp(Rdevs[t-1]-0.5*sigR^2)
        } else if (a>1 && a<nages) {
          Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*Ht[t-1])
        } else if (a==nages) {
          Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*Ht[t-1])+Nat[a, t-1]*exp(-M)*(1.0-Selex_a[a]*Ht[t-1])
        }
      }
      
    }
    
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ##@@ last year derived quantities @@
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SSBt[ntimes+1]=sum(FemaleProp*Nat[,ntimes+1]*Maturity_a*Weight_a)
    availBt[ntimes+1]=sum(Nat[,ntimes+1]*Weight_a*Selex_a)
    Bt[ntimes+1]=sum(Nat[,ntimes+1]*Weight_a)
    
    BStatus=availBt/availB0
    SSBStatus=SSBt/SSB0
    
    It=q*availBt[1:ntimes]*exp(rnorm(ntimes, 0, sd=sqrt(tau2)))
    
    
    if(Bt[ntimes+1]>0 & all(Ht>=0) & all(Ht<=1) ) {break}
    
  }
  
  if(ObtainRef==1) {
    
    bioRef=NumericMSY(nages=nages,
                      R0=R0,
                      h=h,
                      eqSurv=eqSurv,
                      Weight_a=Weight_a,
                      Maturity_a=Maturity_a,
                      FemaleProp=FemaleProp,
                      Selex_a=Selex_a,
                      M=M) 
    
    MSY=bioRef$MSY
    Hmsy=bioRef$Hmsy
    Bmsy=bioRef$Bmsy
    
    BtoverBmsy=availBt/Bmsy
    HtoverHmsy=Ht/Hmsy
    
  } else {
    MSY=NULL
    Hmsy=NULL
    Bmsy=NULL
    BtoverBmsy=NULL
    HtoverHmsy=NULL
  }
  
  
  
  
  
  return(list("Nat"=Nat,
              "availBt"=availBt,
              "SSBt"=SSBt,
              "Bt"=Bt,
              "Ht"=Ht,
              "BStatus"=BStatus,
              "SSBStatus"=SSBStatus,
              "BtoverBmsy"= BtoverBmsy,
              "HtoverHmsy"= HtoverHmsy,
              "M"=M,
              "R0"=R0,
              "h"=h,
              "tau2"=tau2,
              "q"=q,
              "It"=It,
              "Ht"=Ht,
              "MSY"=MSY,
              "Hmsy"=Hmsy,
              "Bmsy"=Bmsy))
}



###################################
######### simulate data ###########
###################################

AspmSimData=function(a_mat, 
                     nages, 
                     sig_mat, 
                     kappa, 
                     Linf, 
                     a0, 
                     lwa, 
                     lwb,
                     a_sel, 
                     sig_sel, 
                     sigR,
                     Ct,
                     ObtainRef=T,
                     ParsPrior) {
  
  ntimes=length(Ct)
  
  
  ## sub-models
  Maturity_a=maturity_at_age(a_mat=a_mat, nages=nages, sig_mat=sig_mat)
  Length_a=length_at_age(kappa=kappa, Linf=Linf, a0=a0, nages=nages)
  Weight_a=weight_at_age(lwa=lwa, lwb=lwb, length_a = Length_a)
  Selex_a=selectivity_at_age(a_sel=a_sel, nages=nages, sig_sel=sig_sel)
  
  
  DynSim=ASPM_dynamics(Maturity_a = Maturity_a,
                       Weight_a = Weight_a,
                       Selex_a = Selex_a,
                       nages=nages,
                       sigR=sigR,
                       ObtainRef=ObtainRef,
                       Ct=Ct,
                       ParsPrior=ParsPrior)
  
  return(DynSim)
  
}

