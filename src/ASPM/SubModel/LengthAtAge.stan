
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