

// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————


//@@@@@@@@@@@@@@@@
//@@@ Maturity @@@
//@@@@@@@@@@@@@@@@

vector Mat(real a_mat, 
           real sig_mat,
           int nages
           ) {
             
    vector[nages] Maturity;
    
    for (a in 1:nages) {
    
    if(a<a_mat) {
    Maturity[a]=0;
    } else if (a>=a_mat) { 
    Maturity[a]=1;
    }
    
    //Maturity[a]=1. ./ (1.+exp( (-(a-a_mat)/sig_mat) ));
    }
    
    
    
    return Maturity;         
  }
  