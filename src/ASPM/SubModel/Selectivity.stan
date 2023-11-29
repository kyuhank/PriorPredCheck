
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
