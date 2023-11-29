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
