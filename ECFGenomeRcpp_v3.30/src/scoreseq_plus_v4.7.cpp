#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//this function assumes that pwm35 and pwm10 have the same lengths
//version 2 adds the functionality of scoreseq_cpp11 to the function!
//version 3 does the two pwms in parallel
//version 4 adds the ability to do several pwms in parallel.
//it expects them as rbound matrices
//version 4.1 fixes a few small mistakes
//version 4.2 introduces combinatoric shuffling in scoreseq_plus_fwd
//version 4.3 is identical to verson 4.2 but has better annotation
//version 4.4 implements manual loop unrolling for number of ecfs 1:5
//Version 4.5 attempts to implement the ListOf<T>
//version 4.6 changes the location of a for loop for a drastic speedup!
//version 4.61 addresses an issue with the default case
//version 4.7 implements list return type from scoreseq_cpp3 (later, different size PWMs)
//and some speedups including loop unrolling for 5 pwms


//this is a function to which R wont or shouldnt have access to but will be used inside this function
//it takes as an input stacked pwm1 and stacked pwm2 stacked here means 4*n of rows, m number of cols
//it then outputs a matrix with cols as the scores of the two pwms in this way 
//pwm1.0, pwm2.0, pwm1.1, pwm2.1, pwm1.2...etc
//as written here, pwm1 and pwm2 have to have the same size
List scoreseq_cpp3(NumericMatrix pwm1, NumericMatrix pwm2, CharacterVector seqs, int number_of_ecfs) {
  
  //get the size of the pwm1 and pwm2
  //for now we will continiue to assume the two pwms are the same size
  int ncol = pwm1.ncol();

  //get the lengt of the sequence
  int seq_len = seqs.size();

  //initialize the output. it will look like
  //there will be two results that are returned as a list
  NumericMatrix final_results_1(seq_len,number_of_ecfs);
  NumericMatrix final_results_2(seq_len,number_of_ecfs);
  
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
            final_results_1((j-i), b) += pwm1((0+(b*4)),i);
            final_results_2((j-i), b) += pwm2((0+(b*4)),i);
          }
        }
        break;
      
      case 'c':
        //for every sequence downstream
        for (int i = 0; i < i_index; i++ ){
          //for every ecf...
          for (int b = 0; b < number_of_ecfs; b++){
            final_results_1((j-i), b) += pwm1((1+(b*4)),i);
            final_results_2((j-i), b) += pwm2((1+(b*4)),i);
          }
        }
        break;
      
      case 'g':
        //for every sequence downstream
        for (int i = 0; i < i_index; i++ ){
          //for every ecf...
          for (int b = 0; b < number_of_ecfs; b++){
            final_results_1((j-i), b) += pwm1((2+(b*4)),i);
            final_results_2((j-i), b) += pwm2((2+(b*4)),i);
          }
        }
        break;
      
      case 't':
        //for every sequence downstream
        for (int i = 0; i < i_index; i++ ){
          //for every ecf...
          for (int b = 0; b < number_of_ecfs; b++){
            final_results_1((j-i), b) += pwm1((3+(b*4)),i);
            final_results_2((j-i), b) += pwm2((3+(b*4)),i);
          }
        }
        break;
      
      default:
        //for every sequence downstream
        for (int i = 0; i < i_index; i++ ){
          //for every ecf, do the accumulate. if the letter is not acgt, then take whatever the minimum
          for (int b = 0; b < number_of_ecfs; b++){
            int b4 = b*4;
            final_results_1((j-i), b) += std::min(std::min(std::min(pwm1((3+b4),i),pwm1((2+b4),i)),pwm1((1+b4),i)),pwm1((0+b4),i));
            final_results_2((j-i), b) += std::min(std::min(std::min(pwm2((3+b4),i),pwm2((2+b4),i)),pwm2((1+b4),i)),pwm2((0+b4),i));
          }
        }
      break;
      }
  }
  
  //return statement
  return Rcpp::List::create(final_results_1(seq(0,(seq_len-ncol)),_),
                            final_results_2(seq(0,(seq_len-ncol)),_));
  
}

//this function combines the pwm1 and pwm2 score here refered to as pwm35 and pwm10
//in all possible combinations with a spacer
//the spacer is the distance from the start of one motif to the start of the second
// [[Rcpp::export]]
List scoreseq_plus_combinatoric2(NumericMatrix pwm35, NumericMatrix pwm10, ListOf<CharacterVector> sequences, IntegerVector spacer) {

  //howmany stacked motifs are there?
  int number_of_ecfs = pwm35.nrow()/4;
  
  //get the number of sequences and 
  int number_of_sequences = sequences.size();
  
  //set up the results vectors
  //there is an all_scores, all_dist, and all_spacers
  NumericMatrix all_scores(number_of_sequences,(number_of_ecfs*number_of_ecfs));
  NumericMatrix all_dist(number_of_sequences,(number_of_ecfs*number_of_ecfs));
  std::fill(all_dist.begin(), all_dist.end(), 0);
  NumericMatrix all_spacers(number_of_sequences,(number_of_ecfs*number_of_ecfs));
  std::fill(all_spacers.begin(), all_spacers.end(), spacer[0]);
    
  for (int a = 0; a < number_of_sequences; a++){
  
  //first get the pwm35 scores and pwm10 scores
  List pwm_scores = scoreseq_cpp3(pwm35, pwm10, sequences[a], number_of_ecfs);
  NumericMatrix pwm_scores_1 = pwm_scores[0];
  NumericMatrix pwm_scores_2 = pwm_scores[1];
  
  //initialize some basic parameters
  int pwm_scores_size = pwm_scores_1.nrow();

  //for each shuffle of -10 and -35
  for (int pwm35_pos = 0; pwm35_pos < number_of_ecfs; pwm35_pos++){
    for (int pwm10_pos = 0; pwm10_pos < number_of_ecfs; pwm10_pos++){    
    
    //this is where in the results file to put them - a non duplicating index
    int index_of_results = (pwm35_pos*number_of_ecfs)+pwm10_pos;
    
    //initialize these values 
    all_scores(a,index_of_results) = pwm_scores_1(0,pwm35_pos) + pwm_scores_2(spacer[0],pwm10_pos);
    
    //we can do this as a nested for loop!
    for (const auto& spac : spacer){   
    
      //length of the second loop depends on the spacer we are using (shorter spacer means slightly more values)
      int sum_scores_length = (pwm_scores_size-spac);
      
      //for each thing in each spacer length
      for (int i = 0; i < sum_scores_length; i++){
        
        //current value is caluclated as the sum of the appropriate pwm_scores
        double current_value = pwm_scores_1(i,pwm35_pos) + pwm_scores_2((i+spac),pwm10_pos);
        
        //if the current value is higher than the previous highest value
        if (current_value > all_scores(a,index_of_results)){
          
          //replace all of the info about the previous value with all of the info about the new value
          all_scores(a,index_of_results) = current_value;
          all_dist(a,index_of_results) = i;
          all_spacers(a, index_of_results) = spac;
          
        }    
      }   
        
    }
  }
  }
  }
  
  //convert all_dist to R indeces
  all_dist = all_dist+1;
  
  //return a list of things
  return Rcpp::List::create(Rcpp::Named("scores") = all_scores,
                            Rcpp::Named("dist") = all_dist,
                            Rcpp::Named("spacers") = all_spacers);
  
}