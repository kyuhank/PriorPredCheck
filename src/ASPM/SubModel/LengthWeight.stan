
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
