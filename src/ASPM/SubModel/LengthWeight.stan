

// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————


//@@@@@@@@@@@@@@@@@@@@@
//@@@ Length-Weight @@@
//@@@@@@@@@@@@@@@@@@@@@

vector Weight(vector Length_a,
              real lwa, 
              real lwb,
              int nages) {
    
    vector[nages] LW=lwa*Length_a^(lwb);

    return LW;         
  }
