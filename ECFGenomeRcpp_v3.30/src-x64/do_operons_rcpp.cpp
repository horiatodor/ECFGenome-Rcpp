#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix do_operon_Rcpp(NumericVector all_scores, List operon_two, IntegerVector matches_to_all_scores) {

  //get the number of potential operons
  int number_of_operons = operon_two.length();

  //results 
  NumericVector temp_results(clone(matches_to_all_scores));
  NumericVector safe_matches(clone(matches_to_all_scores));
  
  //for each potential operon
  for (int i = 0; i < number_of_operons; i++){
  
    //get the hits
    NumericVector ops_of_int = operon_two[i];

    //get the number of genes in the potential operon
    int number_of_genes_in_operon = ops_of_int.length();
    
    //and current top scores
    double current_highest_score = all_scores((matches_to_all_scores(ops_of_int(0))-1));
    int current_highest_index = matches_to_all_scores(ops_of_int(0));
    
    //loop starting with 1 (equivalent to 2)
    for (int j = 1; j < number_of_genes_in_operon; j++){
      
      //check the score, if its bigger...
      if (current_highest_score > all_scores((matches_to_all_scores(ops_of_int(j))-1))){
        temp_results(ops_of_int(j)) = current_highest_index;
      } else {
        current_highest_score = all_scores((matches_to_all_scores(ops_of_int(j))-1));
        current_highest_index = matches_to_all_scores(ops_of_int(j));
      }
      
    }
    
  }
  
  //
  LogicalVector to_keep = safe_matches != temp_results;
    
  safe_matches = safe_matches[to_keep];
  temp_results = temp_results[to_keep];
  
  NumericMatrix results(safe_matches.length(),2);
  
  results(_,0) = safe_matches;  
  results(_,1) = temp_results;  
  
  //
  return results;
  
}
