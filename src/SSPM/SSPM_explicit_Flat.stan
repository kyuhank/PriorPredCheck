
// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————


data {
  int ntimes;
  vector[ntimes] It;
  vector[ntimes] Ct;
  
  vector[2] PriorMu_logKr;
  matrix[2, 2] PriorSigma_logKr;
  vector[2] sigma2Prior;
  vector[2] tau2Prior;
  vector[2] qPrior;
  matrix[2, ntimes] HPrior;
  }
  
transformed data {
  //@@ cholesky decomposition
  matrix[2,2] L=cholesky_decompose(PriorSigma_logKr);
}

parameters {
  vector[2] KrRaw;
  real<lower=log(qPrior[1]), upper=log(qPrior[2]) > logq;
  real<lower=0> sigma2;
  real<lower=0> tau2;
  
  vector<lower=0>[ntimes+1] Pt;
  vector<lower=0, upper=1>[ntimes] Ht;
}

transformed parameters {
 
 vector[2] KrPars=PriorMu_logKr+L*KrRaw;
 
 real K=exp(KrPars[1]);
 real r=exp(KrPars[2]);
 real MSY=(r*K)/4;
 real Hmsy=r/2;
 real Bmsy=K/2;
 
 vector[ntimes] HtoverHmsy;
 vector[ntimes+1] BtoverBmsy;
 vector[ntimes+1] BStatus;
 
 vector<lower=0>[ntimes+1] PtMedian;
 vector[ntimes+1] Bt;
 real q=exp(logq);
 vector<lower=0>[ntimes] PredCt;
 vector[ntimes] Cdiff;
 
 PtMedian[1]=1;
 Bt[1]=Pt[1]*K;
 
 for(t in 2:(ntimes+1)) {
   
   PredCt[t-1]=Bt[t-1]*Ht[t-1];
   
   PtMedian[t]=Pt[t-1]+r*(1.0-Pt[t-1])*Pt[t-1]-PredCt[t-1]/K;
   Bt[t]=Pt[t]*K;
 
 }
 
 HtoverHmsy=Ht/Hmsy;
 BtoverBmsy=Bt/Bmsy;
 
 Cdiff=(log(Ct)-log(PredCt))/0.01;
 
 BStatus=Bt/K;
 
}



model {
  
    KrRaw~ std_normal();
  
    // Priors specified in Meyer and Millar 1999  
    sigma2 ~ inv_gamma(sigma2Prior[1], sigma2Prior[2]);
    tau2 ~ inv_gamma(tau2Prior[1], tau2Prior[2]);
    
    // Process likelihood
    Pt ~ lognormal(log(PtMedian), sqrt(sigma2));
    
    // Observation likelihood
   for(t in 1:ntimes) {
      target+=normal_lpdf( log(It[t]) | logq +log(K)+log(Pt[t]), sqrt(tau2));
    }
    
    Cdiff~std_normal();
    
  for(t in 1:ntimes) {
    Ht[t]~beta(HPrior[1,t],HPrior[2,t]);
  }
}


generated quantities {
  
 vector[2] KrRawPPC;
 vector[2] KrPPC;
 real sigma2PPC;
 real KPPC;
 real rPPC;
 
 vector[ntimes+1] proError;
 vector[ntimes+1] PtPPC;
 vector[ntimes+1] BtPPC;
 
 vector[ntimes] HtPPC;
 vector[ntimes] CtPredPPC;
 
 for(i in 1:2) {
   KrRawPPC[i]=normal_rng(0, 1);
 }
 
 KrPPC=PriorMu_logKr+L*KrRawPPC;
 
 KPPC=exp(KrPPC[1]); 
 rPPC=exp(KrPPC[2]);
 
 for(i in 1:ntimes) { 
 HtPPC[i]=beta_rng(HPrior[1,i],HPrior[2,i]);
 }
 
 sigma2PPC=inv_gamma_rng(sigma2Prior[1], sigma2Prior[2]);
 proError[1]=normal_rng(0, sqrt(sigma2PPC));

 PtPPC[1]=exp(proError[1]);
 BtPPC[1]=PtPPC[1]*KPPC;

 for(t in 2:(ntimes+1) ) {
  CtPredPPC[t-1]=BtPPC[t-1]*HtPPC[t-1];
  proError[t]=normal_rng(0, sqrt(sigma2));
  PtPPC[t]=exp(log(PtPPC[t-1]+rPPC*(1.0-PtPPC[t-1])*PtPPC[t-1]-CtPredPPC[t-1]/KPPC)+proError[t]);
  BtPPC[t]=PtPPC[t]*KPPC;
 }
  
}

