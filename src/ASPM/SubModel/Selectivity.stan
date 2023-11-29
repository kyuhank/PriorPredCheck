
// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————



//@@@@@@@@@@@@@@@@@@@
//@@@ Selectivity @@@
//@@@@@@@@@@@@@@@@@@@
  
vector Selex(real a_sel, 
             real sig_sel,
             int nages
             ) {
    
    vector[nages] Selec;
    
    for(a in 1:nages) {
    Selec[a]=1. ./ (1.+exp( (-(a-a_sel)/sig_sel) ));
    }
    
    return Selec;         
  }
