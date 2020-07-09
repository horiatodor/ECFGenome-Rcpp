#include <Rcpp.h>
using namespace Rcpp;
//this function does the thing for a single pwm
//v2 can do multiple at a time to take advantage of effciencies 
//v2.1 uses a list construct for sequences
//v3 uses a much better sequence walking algorithm

//this is a function to which R wont or shouldnt have access to but will be used inside this function
//it takes as an input stacked pwm1 and stacked pwm2 stacked here means 4*n of rows, m number of cols
//it then outputs a matrix with cols as the scores of the two pwms in this way 
//pwm1.0, pwm1.1, pwm1.2...etc
NumericMatrix scoreseq_cpp1(NumericMatrix pwm1, CharacterVector seqs, int number_of_ecfs) {
  
  //get the size of the pwm1.
  int ncol = pwm1.ncol();
  
  //get the length of the sequence
  int seq_len = seqs.size();

  //initialize the output. it will look like
  //pwm1.0, pwm2.0, pwm1.1, pwm2.1, pwm1.2...etc
  NumericMatrix final_results(seq_len,number_of_ecfs);
  
  //for all of the relevant positions in the sequence
  for(int j = 0; j < seq_len ; j++) {
    
    //get the character of the sequence at the relevant position (here, where in the pwm + where in the sequence)
    char curr = seqs[j][0];
    int i_index = std::min(ncol, (j+1));
    
    //now open a switch
    switch(curr) {
    case 'a':
      //for every sequence downstream
      for (int i = 0; i < i_index; i++ ){
        //for every ecf...
        for (int b = 0; b < number_of_ecfs; b++){
          final_results((j-i), b) += pwm1((0+(b*4)),i);
        }
      }
      break;
    case 'c':
      //for every sequence downstream
      for (int i = 0; i < i_index; i++ ){
        //for every ecf...
        for (int b = 0; b < number_of_ecfs; b++){
          final_results((j-i), b) += pwm1((1+(b*4)),i);
        }
      }
      break;
    case 'g':
      //for every sequence downstream
      for (int i = 0; i < i_index; i++ ){
        //for every ecf...
        for (int b = 0; b < number_of_ecfs; b++){
          final_results((j-i), b) += pwm1((2+(b*4)),i);
        }
      }
      break;
    case 't':
      //for every sequence downstream
      for (int i = 0; i < i_index; i++ ){
        //for every ecf...
        for (int b = 0; b < number_of_ecfs; b++){
          final_results((j-i), b) += pwm1((3+(b*4)),i);
        }
      }
      break;
    default:
      //for every sequence downstream
      for (int i = 0; i < i_index; i++ ){
        //for every ecf...
        for (int b = 0; b < number_of_ecfs; b++){
          int b4 = b * 4;
          final_results((j-i), b) += std::min(std::min(std::min(pwm1((3+b4),i),pwm1((2+b4),i)),pwm1((1+b4),i)),pwm1((0+b4),i));
        }
      }
    break;
    }
  }

  //return statement
  return final_results(seq(0,(seq_len-ncol)),_);
  
}

//this function gets some basic info about where in the sequence the best hit is

// [[Rcpp::export]]
List scoreseq_plus_single2(NumericMatrix pwm35, ListOf<CharacterVector> sequences) {

  //howmany stacked motifs are there?
  int number_of_ecfs = pwm35.nrow()/4;
  
  //get the number of sequences and 
  int number_of_sequences = sequences.size();
  
  //set up the results vector
  NumericMatrix all_scores(number_of_sequences,number_of_ecfs);
  NumericMatrix all_dist(number_of_sequences,number_of_ecfs);
  NumericMatrix all_spacers(number_of_sequences,number_of_ecfs);
  
  for (int a = 0; a < number_of_sequences; a++){
    
  //first get the pwm35_scores
  NumericMatrix pwm_scores = scoreseq_cpp1(pwm35, sequences[a], number_of_ecfs);
    
  //initialize some basic parameters
  int pwm_scores_size = pwm_scores.nrow();

  //for each shuffle of -10 and -35
  for (int index_of_results = 0; index_of_results < number_of_ecfs; index_of_results++){
    
    //initialize these values 
    all_scores(a,index_of_results) = pwm_scores(0,index_of_results);
    all_dist(a,index_of_results) = 1;

    //for each thing in each spacer length
    for (int i = 0; i < pwm_scores_size; i++){
    
      //current value is caluclated as the sum of the appropriate pwm_scores
      //double current_value = pwm_scores(i,index_of_results);
      
      //if the current value is higher than the previous highest value
      if (pwm_scores(i,index_of_results) > all_scores(a,index_of_results)){
        
        //replace all of the info about the previous value with all of the info about the new value
        all_scores(a,index_of_results) = pwm_scores(i,index_of_results);
        all_dist(a,index_of_results) = (i+1);

      }    
    }  
  }
  }
  
  //return the matrix (description is above)
  return Rcpp::List::create(Rcpp::Named("scores") = all_scores,
                            Rcpp::Named("dist") = all_dist);

}