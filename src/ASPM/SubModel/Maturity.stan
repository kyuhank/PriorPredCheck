
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
  