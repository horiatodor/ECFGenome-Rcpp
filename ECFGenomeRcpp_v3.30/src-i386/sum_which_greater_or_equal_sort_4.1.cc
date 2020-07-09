// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//this function sorts but otherwise performs identically to sum_which_greater
//its ~2-3x faster!
//v4 renames the variables, also sorts z
//of_int = x, to_sum = y, avector = z

//v4.1 must have x in decreasing order to work, but is faster

// [[Rcpp::depends("RcppArmadillo")]]
Rcpp::NumericVector ordera(arma::vec x) {
  return(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(arma::sort_index( x ))) );
}

// [[Rcpp::export]]
NumericVector sum_which_greater_or_equal(NumericVector x, NumericVector y, NumericVector z) {
  
  //initialize the return vector
  int z_size = z.size();
  int x_size = x.size();
  NumericVector results(x_size);
  NumericVector orderz = rev(ordera(z));
  
  //the loops
  //for everything in z
  int start_here = 0;
  
  for(int i = 0; i < z_size; i++) {
    
    double comp_thing = z[orderz[i]];
    
    for(int j = start_here; j < x_size; j++) {
      
      if (x[j] <= comp_thing){
        
        results[Range(j,(x_size-1))] = results[Range(j,(x_size-1))] + y[orderz[i]]; 
        
        start_here = j;
        
        break;
        
      }
      
    }
    
  }
  
  // return the thing
  return results;
  
}