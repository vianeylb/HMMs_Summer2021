// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <omp.h>

static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = false,
                         int const cores = 1) {  
  using arma::uword;
  omp_set_num_threads(cores);
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
#pragma omp parallel for schedule(static) private(z)
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);   
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}
=======
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <omp.h>

static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = false,
                         int const cores = 1) {  
  using arma::uword;
  omp_set_num_threads(cores);
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat R;
  bool success;
  success = arma::chol(R, sigma);
  if (success == false)
  {
    return out;
  }
  else{
    arma::mat const rooti = arma::inv(trimatu(R));
    double const rootisum = arma::sum(log(rooti.diag())), 
      constants = -(double)xdim/2.0 * log2pi, 
      other_terms = rootisum + constants;
    
    arma::rowvec z;
#pragma omp parallel for schedule(static) private(z)
    for (uword i = 0; i < n; i++) {
      z = (x.row(i) - mean);
      inplace_tri_mat_mult(z, rooti);   
      out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
    
    if (logd)
      return out;
    return exp(out);
  }
  
}
  