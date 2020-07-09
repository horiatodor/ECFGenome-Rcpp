// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sum_which_greater(NumericVector of_int, NumericVector to_sum, NumericVector avector) {
  
  //initialize the return vector
  int avector_size = avector.size();
  int of_int_size = of_int.size();
  NumericVector results(of_int_size);
  
  //the loops
  for(int i = 0; i < avector_size; i++) {
    
    double comp_thing = avector[i];
    
    for(int j = 0; j < of_int_size; j++) {
      
      if (of_int[j] < comp_thing){
        results[j] += to_sum[i];
      }
      
    }
    
  }
  
  // return the thing
  return results;
  
}