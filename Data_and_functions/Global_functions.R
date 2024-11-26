
# ———————————————————————————————
# Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# ———————————————————————————————


GetMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


GetPval=function(StanFitObj, QuantOfInterest, TrueVal) {
  
  Draws=StanFitObj$draws(QuantOfInterest)
  
  nSamps=dim(Draws)[1]*dim(Draws)[2]
  
  apply(sweep(as_draws_matrix(Draws), 2, TrueVal, ">"), 2, sum)/nSamps
  
}

bounded_inv_logit=function(y, Lower, Upper) {
  
  return((Upper-Lower)/(1.0+exp(-y))+Lower)
  
}


bounded_logit=function(x, Lower, Upper) {
  
  return(-log(((Upper-Lower)/(x-Lower))-1))
  
}


logit=function(p) {
  return(log(p/(1 - p)) )
}


lnSigmaFromCV=function(CV) {
  sqrt(log(CV^2 + 1))
}




SBC=function(StanObj,
             QuantOfInterest,
             ModelType,
             InputData, 
             InitPars,
             seed=123,
             nchains=5,
             nwarmup=10000,
             nsample=20000,
             adaptDelta=0.99,
             maxTree=10,
             nsim=200,
             nPost=100,
             verbose=F,
             nrefresh=0,
             TailMinEffSamp=100,
             BulkMinEffSamp=100,
             ParsPrior,
             ...) {
  
  trueVals=list()
  FitSum=list()
  FitSumPML=list()
  diagnostic=list()
  Pvals=list()
  postSamples=list()
  rankStat=list()
  
  TailEffSamples=list()
  BulkEffSamples=list()
  
  Ct=InputData$Ct
  
  Cdiff=c()
  
  
  if(ModelType=="ASPM") {
    a_mat=InputData$a_mat 
    nages=InputData$nages 
    sig_mat=InputData$sig_mat 
    kappa=InputData$kappa 
    Linf=InputData$Linf 
    a0=InputData$a0 
    lwa=InputData$lwa
    lwb=InputData$lwb 
    a_sel=InputData$a_sel
    sig_sel=InputData$sig_sel
    sigR=InputData$sigR 
  }
  
  
  for ( i in 1:nsim) {
    
    print(paste0("iter-",i))
    
    if(ModelType=="SSPM") {
      SimData=SspmSimData(ParsPrior,
                          ntimes=ntimes,
                          Ct=Ct)
    }
    
    if (ModelType=="ASPM") {
      SimData=AspmSimData(a_mat=a_mat, 
                          nages=nages, 
                          sig_mat=sig_mat, 
                          kappa=kappa, 
                          Linf=Linf, 
                          a0=a0, 
                          lwa=lwa, 
                          lwb=lwb, 
                          a_sel=a_sel, 
                          sig_sel=sig_sel,
                          sigR = sigR, 
                          Ct = Ct,
                          ParsPrior=ParsPrior) 
    }
    
    InputData$It=SimData$It
    
    
    trueVals[[i]]=unlist(SimData[c(QuantOfInterest, "MSY", "Hmsy", "Bmsy", "BtoverBmsy", "HtoverHmsy")])
    
    capture.output({SimFit=StanObj$sample(data=InputData, 
                                          init = rep(InitPars, nchains),
                                          seed = seed, 
                                          chains =  nchains, 
                                          parallel_chains =  nchains,  
                                          iter_warmup = nwarmup,
                                          iter_sampling = nsample,
                                          adapt_delta = adaptDelta,
                                          max_treedepth = maxTree, 
                                          refresh = nrefresh,
                                          show_messages = verbose)}, file=nullfile())
    
    
    if (ModelType=="ASPM") {
      
      capture.output({GenQuant=StanObj$generate_quantities(fitted_params = SimFit, data=InputData, parallel_chains = nchains)}, file=nullfile())
      
     DrawsSamples=cbind(as_draws_matrix(SimFit$draws( QuantOfInterest )), as_draws_matrix(GenQuant$draws( c("MSY", "Hmsy", "Bmsy", "BtoverBmsy", "HtoverHmsy") )))
    
    #  DrawsSamples=as_draws_matrix(SimFit$draws( c(QuantOfInterest, "MSY", "Hmsy", "Bmsy", "BtoverBmsy", "HtoverHmsy") ))
      } 
    
    if (ModelType=="SSPM") {
      DrawsSamples=as_draws_matrix(SimFit$draws(c(QuantOfInterest, "MSY", "Hmsy", "Bmsy", "BtoverBmsy", "HtoverHmsy") ))
    } 
    
    
    postSamples[[i]]=DrawsSamples[seq(1,  nsample*nchains, length.out=nPost),]
    
    Cdiff[i]=sum(c(SimFit$draws("Cdiff")))
    
    TailEffSamples[[i]]=apply(DrawsSamples, 2, ess_tail)
    BulkEffSamples[[i]]=apply(DrawsSamples, 2, ess_bulk)
    
    rankStat[[i]]=apply(sweep(postSamples[[i]], 2, trueVals[[i]], FUN="<"), 2, sum)
    
    
    FitSum[[i]]=summary(as_draws_df(DrawsSamples),'mean', ~quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T),'sd',default_convergence_measures(), default_mcse_measures())
    
    
    diaSum=SimFit$diagnostic_summary()  
    diagnostic[[i]]=c("num_divergent"=sum(diaSum$num_divergent), "num_max_treedepth"=sum(diaSum$num_max_treedepth), "rhat"=FitSum[[i]]$rhat)
    
  }
  
  trueVals=bind_rows(trueVals)
  FitSum=bind_rows(FitSum, .id="SimIter")
  
  
  EstMedian=matrix(FitSum$`50%`, nrow=nsim, ncol=dim(trueVals)[2], byrow = T)
  EstMean=matrix(FitSum$mean, nrow=nsim, ncol=dim(trueVals)[2], byrow = T)
  
  TailEffSamples=do.call(rbind, TailEffSamples)
  BulkEffSamples=do.call(rbind, BulkEffSamples)
  
  diagnostic=bind_rows(diagnostic, .id="SimIter")
  
  diagnoPass=unlist(diagnostic[,2])==0 & apply(diagnostic[,-c(1:3)], 1, function(x) all(x<1.01) ) 
  diagnoPass[is.na(diagnoPass)]=T
  
  TailEffPass=apply(TailEffSamples, 1, function(x) all(x>TailMinEffSamp))
  TailEffPass[is.na(TailEffPass)]=T
  
  BulkEffPass=apply(BulkEffSamples, 1, function(x) all(x>BulkMinEffSamp))
  BulkEffPass[is.na(BulkEffPass)]=T
  
  EfftrueVals=bind_rows(trueVals[which(diagnoPass & TailEffPass & BulkEffPass),])
  rankStat=do.call(rbind,rankStat[which(diagnoPass & TailEffPass & BulkEffPass)])
  ExpectedPost=do.call(rbind, postSamples[which(diagnoPass & TailEffPass & BulkEffPass)])
  
  
  nConverged=sum(diagnoPass & TailEffPass & BulkEffPass)
  
  lowB=matrix(FitSum$`2.5%`, nrow=nsim, ncol=dim(trueVals)[2], byrow = T)
  upB=matrix(FitSum$`97.5%`, nrow=nsim, ncol=dim(trueVals)[2], byrow = T)
  
  coverage95=apply(lowB[which(diagnoPass & TailEffPass & BulkEffPass),] < trueVals[which(diagnoPass & TailEffPass & BulkEffPass),] &
                     upB[which(diagnoPass & TailEffPass & BulkEffPass),] > trueVals[which(diagnoPass & TailEffPass & BulkEffPass),], 2, sum)/nConverged
  
  
  EffEstMedian= EstMedian[which(diagnoPass & TailEffPass & BulkEffPass),]
  EffEstMean= EstMean[which(diagnoPass & TailEffPass & BulkEffPass),]
  
  return(list("trueVals"=trueVals,
              "EfftrueVals"=EfftrueVals,
              "FitSum"=FitSum,
              "EstMedian"=EstMedian,
              "EffEstMedian"=EffEstMedian,
              "EstMean"=EstMean,
              "EffEstMean"=EffEstMean,
              "diagnostic"=diagnostic,
              "nConverged"=nConverged,
              "coverage95"=coverage95,
              "rankStat"=rankStat,
              "ExpectedPost"=ExpectedPost,
              "TailEffSamples"=TailEffSamples,
              "BulkEffSamples"=BulkEffSamples,
              "Cdiff"=Cdiff,
              "nSim"=nsim))
  
}
