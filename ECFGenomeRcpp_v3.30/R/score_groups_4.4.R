################################
#this function scores the normal, scores the shuffled 
#

#v4
#this version implements Z-scores 
#instead of pnorm scores

#v4.1
#adds function to get max pwm scores
#fixed issue when no na 

#v4.2
#added input to deal when there is no annotation in score_matrix
#changed calculate distance to only calculate half of the correlations (since the top part of the square
#matrix is identical to the bottom...)

#v4.4 
#cleaned up unused things in score_matrix - specifically fixed the problem with returning NAs
#got rid of the function that does the correlation scoring.

#as currently written NA values are ignored and assigned a value of 1
#to change would also need to multiply
score_matrix <- function(amatrix, weights_of_species, is.annotated = TRUE){
  
  #the thing
  last_col <- ifelse(is.annotated, (dim(amatrix)[2]-1), dim(amatrix)[2])
  matrix_scores <- matrix(NA, dim(amatrix)[1], (last_col-3))

  # were going to go by column
  for (j in 4:last_col){
    
    #all the stuff thats NA 
    na_vals <- which(is.na(amatrix[,j]))
    of_int <- seq(1:dim(amatrix)[1])
      
    #plug in navals, get the other stuff
    if (length(na_vals)!=0){
      matrix_scores[na_vals,(j-3)] <- 0
      scores_of_int <- as.numeric(amatrix[-na_vals,j])
      of_int <- of_int[-na_vals]
    } else {
      scores_of_int <- as.numeric(amatrix[,j])
    }

    #multiply the scores by the appropriate weights
    matrix_scores[of_int,(j-3)] <- weights_of_species[(j-3)]*scores_of_int
    
  }
  
  #return the scores
  return(rowMeans(matrix_scores))
  
}


################################
#get max pwm_scores
get_max_pwm <- function(list_of_pwm_35, list_of_pwm_10, parameters, species_weights){
  
  #
  results <- rep(NA, length=length(list_of_pwm_10))
  
  #for each ECF
  for (i in 1:length(results)){
    
    results[i] <- sum(apply(list_of_pwm_35[[i]],2,max)) + sum(apply(list_of_pwm_10[[i]],2,max)) 
    
  }
  
  #
  return(mean(((results-parameters[1,])/parameters[2,])*species_weights))
  
}