data {
  int ntimes;
  vector[ntimes] It;
  vector[ntimes] Ct;
  
  vector[2] PriorMu_logKr;
  matrix[2, 2] PriorSigma_logKr;
  vector[2] sigma2Prior;
  vector[2] tau2Prior;
  vector[2] qPrior;
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
 vector[ntimes] Ht;
 real q=exp(logq);
 real<lower=0> PtMedTemp;
 vector[ntimes] PredCt;
 vector[ntimes] Cdiff;
 
 PtMedian[1]=1;
 Bt[1]=Pt[1]*K;
 
 for(t in 2:(ntimes+1)) {
   
   
   PtMedTemp=Pt[t-1]+r*(1.0-Pt[t-1])*Pt[t-1];

   if( PtMedTemp<=(Ct[t-1]/K) ) {
     PredCt[t-1]=K*PtMedTemp*0.99;
   } else {
     PredCt[t-1]=Ct[t-1];
   }
   
   PtMedian[t]=Pt[t-1]+r*(1.0-Pt[t-1])*Pt[t-1]-PredCt[t-1]/K;
   Bt[t]=Pt[t]*K;
   
   
   Ht[t-1]=PredCt[t-1]/Bt[t-1];
 }
 
 HtoverHmsy=Ht/Hmsy;
 BtoverBmsy=Bt/Bmsy;
 
  Cdiff=(log(Ct)-log(PredCt))/0.1;


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
      target+=normal_lpdf( log(It[t]) | log(q * K * Pt[t]), sqrt(tau2));
    }
    
  //  Cdiff~exponential(100);
    Cdiff~std_normal();
}



generated quantities{
  vector [ntimes+1] BStatusSim =BStatus;
  }

