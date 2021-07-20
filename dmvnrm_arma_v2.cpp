// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

static double const log2pi = std::log(2.0 * M_PI);
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat const &x,  
                      arma::rowvec const &mean,  
                      arma::mat const &sigma, 
                      bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  // inserting an if statement to address issues with cholesky decomposition
  arma::mat R;
  bool success;
  success = arma::chol(R, sigma);
  if (success == false)
  {
    return out;
  }
  else{
    
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z      = (x.row(i) - mean) * rooti;    
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
  }
}