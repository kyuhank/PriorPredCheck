
// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————



//@@@ include SubModels

functions{
  #include "../SubModel/LengthWeight.stan"
  #include "../SubModel/LengthAtAge.stan"
  #include "../SubModel/Maturity.stan"
  #include "../SubModel/Selectivity.stan"
  #include "../SubModel/StockRecruitment.stan"
  #include "../SubModel/eqQuantities.stan"
  #include "../SubModel/DynamicsLoop.stan"
  #include "../SubModel/Helpers.stan"
  #include "../SubModel/NumericMSY.stan"
}

data {
  
  int ntimes;
  int nages;
  vector[ntimes] It;
  vector[ntimes] Ct;
  
  vector[2] PriorMu_transMh;
  vector[2] PriorSigma_transMh;
  
  matrix[3, 3] Sigma_MR0hRaw;
  
  vector[2] R0Prior;
  vector[2] tau2Prior;
  vector[2] qPrior;
  real sigR;
  real FemaleProp;
  
  // fixed pars for submodels
  real a_mat;
  real sig_mat;
  
  real a_sel;
  real sig_sel;
  
  real Linf;
  real kappa;
  real a0;
  
  real lwa;
  real lwb;
  
}

transformed data {
  //@@ cholesky decomposition
  matrix[3,3] L=cholesky_decompose(Sigma_MR0hRaw);
  
  vector[nages] Maturity_a=Mat(a_mat, sig_mat, nages);
  vector[nages] Selex_a=Selex(a_sel, sig_sel, nages);
  vector[nages] Length_a=Length(Linf, kappa, a0, nages);
  vector[nages] Weight_a=Weight(Length_a, lwa, lwb, nages);
  }

parameters {
  vector[3] MR0hRaw;
  real<lower=log(qPrior[1]), upper=log(qPrior[2])> logq;
  real<lower=0> tau2;
  vector[ntimes] logRt;
}

transformed parameters {
  vector[3] mvnMR0hRaw=L*MR0hRaw;
  
  real M=exp( mvnMR0hRaw[1]*PriorSigma_transMh[1]+PriorMu_transMh[1]);
  real R0=exp( normal_cdf(mvnMR0hRaw[2] | 0, 1)*(log(R0Prior[2])-log(R0Prior[1]))+log(R0Prior[1]) ) ;
  real h=bounded_inv_logit(mvnMR0hRaw[3]*PriorSigma_transMh[2]+PriorMu_transMh[2], 0.2, 1);
  real q=exp(logq);
  
  vector[ntimes] Rt=exp(logRt);
  
  vector[nages] eqSurv=eqSurvRate(M, nages);
  
  matrix[nages+8, ntimes+1] DynOutput= DynamicLoop(nages,
                                                   R0,
                                                   h,
                                                   ntimes,
                                                   eqSurv,
                                                   Weight_a,
                                                   Maturity_a,
                                                   FemaleProp,
                                                   Ct,
                                                   Selex_a,
                                                   M,
                                                   Rt,
                                                   sigR);
                                                
   matrix<lower=0>[nages, ntimes+1] Nat=DynOutput[1:nages,];
   vector<lower=0>[ntimes+1] availBt=DynOutput[nages+1,]';
   vector<lower=0>[ntimes+1] SSBt=DynOutput[nages+2,]';
   vector<lower=0>[ntimes+1] Bt=DynOutput[nages+3,]';
   vector<lower=0>[ntimes] Ht=DynOutput[nages+4,1:ntimes]';
   vector<lower=0>[ntimes] PredCt=DynOutput[nages+5,1:ntimes]';
   vector<lower=0>[ntimes+1] BStatus=DynOutput[nages+6,]';
   vector<lower=0>[ntimes] RtMedian=DynOutput[nages+7,1:ntimes]';
  
   vector<lower=0>[ntimes+1] SSBStatus=DynOutput[nages+8,]';
  
   vector[ntimes] Cdiff=(log(Ct)-log(PredCt))/0.1;
  
  
  vector[ntimes] Rdevs=(logRt-(log(RtMedian)-0.5*sigR^2));
  
  
}

model {
  
  //@@ prior
  //non-centered parameterisation for M, R0, and h
  MR0hRaw~std_normal();
  
  // obs variance for CPUE
  tau2 ~ inv_gamma(tau2Prior[1], tau2Prior[2]);

  // normally distributed logRt
  logRt ~ normal(log(RtMedian)-0.5*sigR^2, sigR);
  
  //@@ likelihood for CPUE data
  for (t in 1:ntimes) {
   target+=normal_lpdf( log(It[t]) | logq+log(availBt[t]), sqrt(tau2));
  }
  
  //@@ penalty catch
   Cdiff~std_normal();
}

generated quantities {
  
  vector[3] refs=NumericMSY(nages,
                            R0,
                            h,
                            eqSurv,
                            Weight_a,
                            Maturity_a,
                            FemaleProp,
                            Selex_a,
                            M);
  
  real MSY=refs[1];
  real Hmsy=refs[2];
  real Bmsy=refs[3];
  
  vector[ntimes+1] BtoverBmsy=availBt/Bmsy;
  vector[ntimes] HtoverHmsy=Ht/Hmsy;
  
}

