
// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————



//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@ Length-at-age (von Bert) @@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

vector Length(real Linf,
                real kappa, 
                real a0,
                int nages) {
    
    vector[nages] LA;
    
    for (a in 1:nages) {
    LA[a]=Linf*(1. - exp(-kappa *(a-a0 )) );
    }

    return LA;
  }