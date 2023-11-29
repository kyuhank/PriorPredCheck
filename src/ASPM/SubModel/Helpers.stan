real bounded_inv_logit(real y, real Lower, real Upper) {
  
  real x=(Upper-Lower)/(1.0+exp(-y))+Lower;
  
  return x;
  
}


real bounded_logit(real x, real Lower, real Upper) {
  
  real y=-log(((Upper-Lower)/(x-Lower))-1);
  
  return y;
  
}


    
vector posfun(real x, real eps) {

real pos_value=0;
real pen=0;

vector[2] value_pen;

if(x>=eps) {
  pos_value=x;
} else {
   pen=0.001*(x-eps)^2.;
   pos_value=eps/(2.-x/eps);           // ADMB
   //pos_value=eps*(1/(1-(x-eps)/eps + (x-eps)^2/eps^2-(x-eps)^3/eps^3 + (x-eps)^4/eps^4 - (x-eps)^5/eps^5   ));  //Ben Bolker's suggestion
}

value_pen[1]=pos_value;
value_pen[2]=pen;

return value_pen; 

}


real normal_copula_lpdf(real u, real v, real rho) {
  real rho_sq = square(rho);
  
  return (0.5 * rho * (-2. * u * v + square(u) * rho + square(v) * rho)) / (-1. + rho_sq)
         - 0.5 * log1m(rho_sq);
}
