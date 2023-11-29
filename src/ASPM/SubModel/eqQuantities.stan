
// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————


//eqSurvivor rate
vector eqSurvRate(real M, 
              int nages) {
              
    vector[nages] SurAtAge=rep_vector(0., nages);
    
    for (a in 1:nages) {
      if(a==1) {
        SurAtAge[a]=1.;
      } else if (a>1 && a<nages) {
        SurAtAge[a]=SurAtAge[a-1]*exp(-M);
      } else if (a==nages) {
        SurAtAge[a]=(SurAtAge[a-1]*exp(-M))/(1. -exp(-M));
      } 
      
    }
    
    return SurAtAge;
    
    }
    