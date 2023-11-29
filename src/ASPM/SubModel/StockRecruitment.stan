
// ———————————————————————————————
// Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// ———————————————————————————————


//@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@ Stock recruitment @@@
//@@@@@@@@@@@@@@@@@@@@@@@@@

//Stock-recruitment relationship
real StockRecruitment(real steepness,
                      real R0,
                      real SSB0,
                      real SSBt) {

      return (4.*R0*steepness*SSBt)/((1.-steepness)*SSB0+(5.*steepness-1.)*SSBt);
      
      }
