
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@ Numerically calculate MSY @@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

vector NumericMSY(int nages,
                  real R0,
                  real h,
                  vector eqSurv,
                  vector Weight_a,
                  vector Maturity_a,
                  real FemaleProp,
                  vector Selex_a,
                  real M
                   ) {
    
    
    int ntimes=200;
    vector[200] eqH;
                
    matrix[nages, ntimes+1] Nat=rep_matrix(0., nages, ntimes+1);
    matrix[nages, ntimes] Cat=rep_matrix(0., nages, ntimes);
    vector[ntimes+1] Bt=rep_vector(0., ntimes+1);
    vector[ntimes+1] availBt=rep_vector(0., ntimes+1);
    vector[ntimes+1] SSBt=rep_vector(0., ntimes+1);
    vector[ntimes+1] BStatus=rep_vector(0., ntimes+1);
    vector[ntimes] Ht=rep_vector(1e-5, ntimes);
    vector[nages] eqNa=R0*eqSurv;
    matrix[ntimes, 200] Ct=rep_matrix(0., ntimes, 200);
    
    real SSB0=0.;
    real B0=0.;
    real availB0=0.;
    
    eqH[1]=0.0;
    for(i in 2:200) {
      eqH[i]=eqH[i-1]+0.005;
    }
  
  //@@@@@@@@@@@@@@@@@@@
  //@@ Eq quantities @@
  //@@@@@@@@@@@@@@@@@@@
  
   for(a in 1:nages) {
     SSB0+=FemaleProp*eqNa[a]*Maturity_a[a]*Weight_a[a];
     B0+=eqNa[a]*Weight_a[a];
     availB0+=eqNa[a]*Weight_a[a]*Selex_a[a];
   }
    
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@ initial N (i.e., t=1) @@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  Nat[,1]=eqNa;
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@ main loop (from t=2)@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@   
  
  for (i in 1:200) {
    for(t in 2:(ntimes+1) ) {


           SSBt[t-1]=sum(FemaleProp* Nat[,t-1] .* Maturity_a .* Weight_a);

        for (a in 1:nages) {
           if (a==1) {
           Nat[a,t]=StockRecruitment(h, R0, SSB0, SSBt[t-1]);
           } else if (a>1 && a<nages) {
           Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*eqH[i]);
           } else if (a==nages) {
           Nat[a,t]=Nat[a-1, t-1]*exp(-M)*(1.0-Selex_a[a-1]*eqH[i])+Nat[a, t-1]*exp(-M)*(1.0-Selex_a[a]*eqH[i]);
           }
           
           
           Cat[a,t-1]=Nat[a,t-1]*Selex_a[a]*eqH[i]*Weight_a[a];
           
           
        }
    
        Ct[t-1,i]=sum(Cat[,t-1]);    
        
    }
    
    
  }
    
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@@ returning quantities @@@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  real MSY;
  real Hmsy;
  real Bmsy;
  int index;
  
  vector[3] refs;
  
  MSY=max(Ct[ntimes,]);
  
  for(i in 1:200) { if(Ct[ntimes,i]==MSY) index=i; } 
  
  Hmsy=eqH[index];
  Bmsy=MSY/Hmsy;
  
  refs[1]=MSY;
  refs[2]=Hmsy;
  refs[3]=Bmsy;
  
  return refs;
  
  }

