

SspmSimData <- function  ( ParsPrior,
                           ntimes,
                           Ct) {
  
  Pt=c()
  Bt=c()
  It=c()
  Ht=c()
  HtoverHmsy=c()
  BtoverBmsy=c()
  
  repeat{
    
    if(is.null(ParsPrior[["FixedQuant"]] )) {  
      logKr=MASS::mvrnorm(1, ParsPrior[["logKr"]][[1]],  ParsPrior[["logKr"]][[2]])
      
      K=exp(logKr[1])
      r=exp(logKr[2])
      
      q=exp(runif(1, log( ParsPrior[["q"]][[1]]), log( ParsPrior[["q"]][[2]])))
      
      sigma2=1/rgamma(1, ParsPrior[["sigma2"]][[1]], ParsPrior[["sigma2"]][[2]])
      tau2=1/rgamma(1, ParsPrior[["tau2"]][[1]], ParsPrior[["tau2"]][[2]])
      
    } else {
      K=ParsPrior[["FixedQuant"]][["K"]]  
      r=ParsPrior[["FixedQuant"]][["r"]]  
      q=ParsPrior[["FixedQuant"]][["q"]]
      sigma2=ParsPrior[["FixedQuant"]][["sigma2"]]  
      tau2=ParsPrior[["FixedQuant"]][["tau2"]]  
    }
    
    
    Ept=rnorm(ntimes+1, 0, sqrt(sigma2))
    Nut=rnorm(ntimes, 0, sqrt(tau2))
    
    
    Pt[1]=exp(Ept[1])
    Bt[1]=Pt[1]*K
    
    Bmsy=K/2
    Hmsy=r/2
    
    
    for(t in 2:(ntimes+1)) {
      Pt[t]=(Pt[t-1]+r*(1.0-Pt[t-1])*Pt[t-1]-Ct[t-1]/K)*exp(Ept[t])
      Bt[t]=Pt[t]*K
      It[t-1]=q*Bt[t-1]*exp(Nut[t-1])
      Ht[t-1]=Ct[t-1]/Bt[t-1]
    }
    
    if(all(Bt>0) & all(Ht>=0) & all(Ht<=1) ) {break}
    
  }
  
  BtoverBmsy=Bt/Bmsy
  HtoverHmsy=Ht/Hmsy
  Bmsy=K/2
  MSY=(r*K)/4
  Hmsy=r/2
  
  BStatus=Bt/K
  
  return(list("K"=K,
              "r"=r,
              "q"=q,
              "sigma2"=sigma2,
              "tau2"=tau2,
              "MSY"=MSY,
              "Bmsy"=Bmsy,
              "Hmsy"=Hmsy,
              "Pt"=Pt,
              "Bt"=Bt,
              "BStatus"=BStatus,
              "It"=It,
              "Ht"=Ht,
              "HtoverHmsy"=HtoverHmsy,
              "BtoverBmsy"=BtoverBmsy))
  }





SspmForDemon <- function  ( ParsPrior,
                            ntimes,
                            Ct) {
  
  Pt=rep(0, ntimes+1)
  Bt=rep(0, ntimes+1)
  It=rep(0, ntimes)
  Ht=rep(0, ntimes)
  
  if(is.null(ParsPrior[["FixedQuant"]] )) {  
    logKr=MASS::mvrnorm(1, ParsPrior[["logKr"]][[1]],  ParsPrior[["logKr"]][[2]])
    
    K=exp(logKr[1])
    r=exp(logKr[2])
    
    q=exp(runif(1, log( ParsPrior[["q"]][[1]]), log( ParsPrior[["q"]][[2]])))
    
    sigma2=1/rgamma(1, ParsPrior[["sigma2"]][[1]], ParsPrior[["sigma2"]][[2]])
    tau2=1/rgamma(1, ParsPrior[["tau2"]][[1]], ParsPrior[["tau2"]][[2]])
    
  } else {
    K=ParsPrior[["FixedQuant"]][["K"]]  
    r=ParsPrior[["FixedQuant"]][["r"]]  
    q=ParsPrior[["FixedQuant"]][["q"]]
    sigma2=ParsPrior[["FixedQuant"]][["sigma2"]]  
    tau2=ParsPrior[["FixedQuant"]][["tau2"]]  
  }
  
  
  Ept=rnorm(ntimes+1, 0, sqrt(sigma2))
  Nut=rnorm(ntimes, 0, sqrt(tau2))
  
  
  Pt[1]=exp(Ept[1])
  Bt[1]=Pt[1]*K
  
  Bmsy=K/2
  Hmsy=r/2
  
  
  for(t in 2:(ntimes+1)) {
    Pt[t]=(Pt[t-1]+r*(1.0-Pt[t-1])*Pt[t-1]-Ct[t-1]/K)*exp(Ept[t])
    Bt[t]=Pt[t]*K
    It[t-1]=q*Bt[t-1]*exp(Nut[t-1])
    Ht[t-1]=Ct[t-1]/Bt[t-1]
    
    if(Bt[t]<0 | Ht[t-1]<0 | Ht[t-1]>1 )  {
      Bt[t]=0
      Pt[t]=0
      break
    }
    
  }
  
  
  #BtoverBmsy=Bt/Bmsy
  #HtoverHmsy=Ht/Hmsy
  Bmsy=K/2
  MSY=(r*K)/4
  Hmsy=r/2
  
  BStatus=Bt/K
  
  return(list("K"=K,
              "r"=r,
              "q"=q,
              "sigma2"=sigma2,
              "tau2"=tau2,
              "MSY"=MSY,
              "Bmsy"=Bmsy,
              "Hmsy"=Hmsy,
              "Pt"=Pt,
              "Bt"=Bt,
              "BStatus"=BStatus,
              "It"=It,
              "Ht"=Ht #,
              #"HtoverHmsy"=HtoverHmsy,
              #"BtoverBmsy"=BtoverBmsy
  ))
}
