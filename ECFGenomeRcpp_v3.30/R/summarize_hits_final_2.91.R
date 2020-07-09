#
#v2.4 fixes major error in 2.3 with respect to how operons were treated!

#v2.41 
#adds distance to ECF start site

#v2.42
#fixes operon identification

#v2.43 
#minor performance increases

#v2.5 major performance improvements

#version 2.521 returns a negative distance when the ECF and the target gene are on opposite strands

#version 2.522 minor speedups

#version 2.6 - completely rewritten to take into account NOG preponderance when calculating Z scores

#version 2.7 changes the definition of outliers to more than 4sd smaller than mean

#version 2.8
#changes the implementation of the cog correction to better correct for edge effects

#version 2.9
#implements MAD and median instead of mean and SD
#adjusted outlier cutoff to 5sd (since SD is now more robust)

#v2.91
#changed qnorm/pnorm to log.p = TRUE to avoid issues with lack of precision around 1 (high z-scores)

###################################################################################
summarize_hits_final_4 <- function(the_whole_thing, gene_starts, ops, 
                                   all_scores, all_dist_from_tss, all_spacers, 
                                   to_match = NULL, ecf_row, correction_factors = NULL){
  
  #ok, first we want to get the mean and standard deviation minus outliers
  hits_mean <- median(all_scores)
  hits_sd <- mad(all_scores)
  
  #get outliers
  #more than 5sd smaller
  outliers <- which(all_scores < (hits_mean-(hits_sd*5)))

  #then do operons
  #if we want to, lets do operons
  if (!is.null(ops)){
    
    if (dim(ops)[1]!=0){
      
      all_scores[ops[,1]] <- all_scores[ops[,2]]
      all_spacers[ops[,1]] <- all_spacers[ops[,2]]
      all_dist_from_tss[ops[,1]] <- -1000
      gene_starts[ops[,1]] <- gene_starts[ops[,2]]
      
    }
  }
  
  #
  if (length(outliers) > 0){all_scores[outliers] <- hits_mean}
  
  #and convert to z scores
  all_scores <- (all_scores-hits_mean)/hits_sd  

  #ok, so we will order the whole thing in order of highest to lowest score, then in order of distance from start site
  order_of_things <- order(all_scores, all_dist_from_tss, decreasing = TRUE)

  #then running match unique should be exactly what we need
  #to implement highest score per COG and if there are multiple than by closest
  #first lets remove the NAs
  if (is.null(to_match)){
    to_match <- names(table(the_whole_thing[,4]))
  }

  #get distance to ECF start site 
  #here we get rid of the gproNOG row in order to now have this info!
  if (!is.na(ecf_row)){
  
    genes_on_same_contig <- which(the_whole_thing[,2] == the_whole_thing[ecf_row,2])
    
    if (length(genes_on_same_contig) != 0){
      
      #the distance
      distance_on_same_contig <- abs(gene_starts[genes_on_same_contig] - gene_starts[ecf_row])
      
      #make the distace negative if theyre not of the same strand
      strand_direction <- ifelse(the_whole_thing[genes_on_same_contig,1] == the_whole_thing[ecf_row,1], 1, -1)
      
      #multiply them
      the_whole_thing[genes_on_same_contig,5] <- distance_on_same_contig*strand_direction
      
    }
  }
  
  #to_return
  to_return <- cbind(gene_starts, 
                     the_whole_thing[,1], 
                     all_scores, 
                     all_spacers, 
                     the_whole_thing[,2], 
                     rep("CDS", length(all_scores)), 
                     all_dist_from_tss,
                     the_whole_thing[,3:9], 
                     stringsAsFactors = FALSE,
                     row.names = NULL)[order_of_things[match(to_match, the_whole_thing[order_of_things,4])],]
  
  #correction factor application
  if (!is.null(correction_factors)){
    to_return[,3] <- qnorm(pnorm(to_return[,3],log.p = TRUE)*correction_factors,log.p = TRUE) 
  }
  
  #return statement
  return(to_return)
  

}
