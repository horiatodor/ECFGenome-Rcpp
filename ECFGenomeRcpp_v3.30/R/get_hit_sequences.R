#v1.0

#function for returning genome sequences
get_hit_sequences <- function(retrieve_sequences, real_hits, fdrs, list_of_truncated_fasta){
  
  #first, lets get the sequences we are interested in
  fdr_range <- retrieve_sequences[[1]]
  rank_range <- retrieve_sequences[[2]]  
  dist_range <- retrieve_sequences[[3]]  
  
  #get the nogs of interest
  within_fdr<- intersect(which(fdrs >= fdr_range[1]),which(fdrs <= fdr_range[2]))
  within_rank <- rank_range[1]:rank_range[2] 
  within_dist <- intersect(which(real_hits[[3]] >= dist_range[1]),which(real_hits[[3]] <= dist_range[2]))  
  
  #
  nogs_of_interest <- intersect(intersect(within_fdr,within_rank),within_dist)
  
  #if there are no nogs of interest, quit and return null
  if (length(nogs_of_interest) == 0){return(NULL)}
  
  #otherwise, time to do some sequence retrieval
  #initialize the outputs
  minus_35_seqs <- NULL
  minus_10_seqs <- NULL
  
  #lets go by genome
  for (j in 1:length(list_of_truncated_fasta)){
    
    #
    how_much_seq <- list_of_truncated_fasta[[j]][[3]]
    
    #for each nog of interest
    for (i in 1:length(nogs_of_interest)){
      
      #if the NOG exists in that genome
      if (!is.na(real_hits[[5]][nogs_of_interest[i],j]) & real_hits[[6]][nogs_of_interest[i],j] != -1000){
        
        #gene name
        gene_name <- real_hits[[5]][nogs_of_interest[i],j]
        spacer <- real_hits[[7]][nogs_of_interest[i],j]

        #get the sequence
        location_in_annot <- which(list_of_truncated_fasta[[j]][[2]][,1] == gene_name)[1]
        sequence_of_int <- list_of_truncated_fasta[[j]][[1]][[location_in_annot]]
        
        #now to subset the sequence...
        #if its on the minus strand rev_comp the sequence
        if (list_of_truncated_fasta[[j]][[2]][location_in_annot,7] == "-"){
          sequence_of_int <- rev(comp(sequence_of_int))
          plus_one_spot <-  real_hits[[6]][nogs_of_interest[i],j] + (how_much_seq[2] - how_much_seq[1])
        } else {
          plus_one_spot <-  real_hits[[6]][nogs_of_interest[i],j] + (how_much_seq[2] - how_much_seq[1]- 1)
        }
        
        #subset
        minus_35_start <- plus_one_spot - spacer - 11 
        minus_35_stop <- plus_one_spot - spacer + 3  
        
        minus_10_start <- plus_one_spot - 11
        minus_10_stop <- plus_one_spot + 3
        
        #adjust the start and stops to take into account the ends of the sequence
        #and add variables for the number on "ns" to add
        start_ns <- 0
        stop_ns <- 0 
        
        #minus35 start
        if (minus_35_start < 1){
          
          start_ns <- 1-minus_35_start
          minus_35_start <- 0
          
        }
        
        #minus10 end
        if (minus_10_stop > length(sequence_of_int)){
          
          stop_ns <- minus_10_stop - length(sequence_of_int)
          minus_10_stop <- length(sequence_of_int)
          
        }
        
        #
        minus_35_seqs <- rbind(minus_35_seqs, c(rep("n", start_ns), sequence_of_int[minus_35_start:minus_35_stop]))
        minus_10_seqs <- rbind(minus_10_seqs, c(sequence_of_int[minus_10_start:minus_10_stop], rep("n", stop_ns)))
        
      }
      
    }
    
  }
  
  return(list(minus_35_seqs, minus_10_seqs))
  
}
