#function to allow parallelization at the iteration level

#v2.1 passes first three and annotation to reduce rbind calls

#v2.11 assumes no annotation on shuffled matrixes

do_an_iteration <- function(list_of_pwm_35, column_shuffle_35, 
                            list_of_pwm_10, column_shuffle_10,
                            list_of_truncated_fasta, list_of_operons, spacers, 
                            species_weights, original_scores, F3andA = NULL){
  
  shuffled_hits <- annotate_group_2(list_of_pwm_35, list_of_pwm_10, 
                                    list_of_truncated_fasta, list_of_operons, spacers, 
                                    shuffled_order_35 = column_shuffle_35, shuffled_order_10 = column_shuffle_10,
                                    cores_to_use = 1, F3andA)
  
  #do the scoring for each resulting permutation and return the results!
  return(Reduce("+",lapply(1:length(shuffled_hits), function (i) score_an_annotated_group(shuffled_hits[[i]], species_weights, original_scores))))                                  

}

#a mini function so i can use an lapply/reduce statement
score_an_annotated_group <- function(annotated_group, species_weights, original_scores){

  shuffled_scores <- score_matrix(annotated_group[[1]], species_weights, is.annotated = FALSE)
  number_better <- sum_which_greater_or_equal(original_scores, annotated_group[[2]], shuffled_scores)
  return(number_better)
  
}  