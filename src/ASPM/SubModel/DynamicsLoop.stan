
//@@@@@@@@@@@@@@@@@@@@@
//@@@ Dynamics loop @@@
//@@@@@@@@@@@@@@@@@@@@@

matrix DynamicLoop(int nages,
                   real R0,
                   real h,
                   int ntimes,
                   vector eqSurv,
                   vector Weight_a,
                   vector Maturity_a,
                   real FemaleProp,
                   vector Ct,
                   vector Selex_a,
                   real M,
                   vector Rt,
                   real sigR
                   ) {
                
    matrix[nages, ntimes+1] Nat=rep_matrix(0., nages, ntimes+1);
    vector[ntimes+1] Bt=rep_vector(0., ntimes+1);
    vector[ntimes+1] availBt=rep_vector(0., ntimes+1);
    vector[ntimes+1] SSBt=rep_vector(0., ntimes+1);
    vector[ntimes+1] BStatus=rep_vector(0., ntimes+1);
    vector[ntimes+1] SSBStatus=rep_vector(0., ntimes+1);
    vector[ntimes] Ht=rep_vector(1e-5, ntimes);
    vector[ntimes] PredCt=rep_vector(1e-5, ntimes);
    
    vector[ntimes] RtMedian=rep_vector(1e-5, ntimes);
    
    vector[nages] eqNa=R0*eqSurv;
    
    real SSB0=0.;
    real B0=0.;
    real availB0=0.;
  
  //@@@@@@@@@@@@@@@@@@@
  //@@ Eq quantities @@
  //@@@@@@@@@@@@@@@@@@@
  
   for(a in 1:nages) {
     SSB0+=FemaleProp*eqNa[a]*Maturity_a[a]*Weight_a[a];
     B0+=eqNa[a]*Weight_a[a];
     availB0+=eqNa[a]*Weight_a[a]*Selex_a[a];
   }
  
  Nat[,1]=eqNa;
     
  //@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@ main loop (from t=2)@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@   
  
    for(t in 2:(ntimes+1) ) {
      
        for(a in 1:nages) {
           SSBt[t-1]+=FemaleProp*Nat[a,t-1]*Maturity_a[a]*Weight_a[a];
           availBt[t-1]+=Nat[a,t-1]*Weight_a[a]*Selex_a[a];
           Bt[t-1]+=Nat[a,t-1]*Weight_a[a];
        }
       
       if( (Ct[t-1]/availBt[t-1])>=1 ) {
         Ht[t-1]=0.99;
         PredCt[t-1]=availBt[t-1]*Ht[t-1];
       } else {
         Ht[t-1]=Ct[t-1]/availBt[t-1];
         PredCt[t-1]=Ct[t-1];
       }
       
      //Ht[t-1]=exp(log(Ct[t-1]/availBt[t-1]));
      //Ht[t-1]=1.-exp(log(1.-Ct[t-1]/availBt[t-1]));
       
       
        for (a in 1:nages) {
           if (a==1) {
           Nat[a,t]=Rt[t-1];
           RtMedian[t-1]=StockRecruitment(h, R0, SSB0, SSBt[t-1]);
           } else if (a>1 && a<nages) {
           Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*Ht[t-1]);
           } else if (a==nages) {
           Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*Ht[t-1])+Nat[a, t-1]*exp(-M)*(1.0-Selex_a[a]*Ht[t-1]);
           }
        }
        
    }
    
    
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@ last year derived quantities @@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    for(a in 1:nages) {
           SSBt[ntimes+1]+=FemaleProp*Nat[a,ntimes+1]*Maturity_a[a]*Weight_a[a];
           availBt[ntimes+1]+=Nat[a,ntimes+1]*Weight_a[a]*Selex_a[a];
           Bt[ntimes+1]+=Nat[a,ntimes+1]*Weight_a[a];
      }
    
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@@ returning quantities @@@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   BStatus=availBt/availB0;
   SSBStatus=SSBt/SSB0;
   
   matrix[nages+8, ntimes+1] OutPut; 
    
   OutPut[1:nages,]=Nat;
   OutPut[nages+1,]=availBt';
   OutPut[nages+2,]=SSBt';
   OutPut[nages+3,]=Bt';
   OutPut[nages+4,1:ntimes]=Ht';
   OutPut[nages+5,1:ntimes]=PredCt';
   OutPut[nages+6,]=BStatus';
   OutPut[nages+7,1:ntimes]=RtMedian';
   OutPut[nages+8,]=SSBStatus';
    
                
  return OutPut;
  
  }

