#version 5 changes
#these are big big changes, so changed function name to annotate_group_2
#we expect this to do all of the things starting with pwms

#v5.1
#changed to score_groups_2

#v5.2
#chaged annotate hits in the shuffle bit, so that they dont return operons 
#we want the significance of individual hits, not of the whole thing.

#v5.3
#his version uses scan and annotate to do the thing

#v5.4
#this version applies iterative scoring 
#what does this look like?
#1. get 5 SHUFFLED HITS
#2. use these shuffled hits to generate PARAMETERS
#3. use these parameters to generate PERCENTILE SCORES for the shuffled and the original hits
#4. use the percentile scores to generate ORIGINAL SCORES and SHUFFLED SCORES
#5. use the shuffled scores to generate un-corrected P-VALUES for the original scores
#6. set an initial score threshold. the threshold is the value of the highest score 
#   that if the corrected pvalue > 1 under best assumptions for the next set
#7. annotate and score another group using the parameters from #2.   
#8. update the scores and pvalues for things above threshold
#9. update the threshold and iterate

#v5.41
#uses custom FDR function for significance assesment

#v5.5
#uses Z-scores for scoring
#as implemented in score_groups_4

#v5.51
#returns results with additional information
#as implemented in score_groups_4.1

#v5.52
#implements reversed pwm approximately half the time

#v5.6
#implements opeoron removal for fdr calculation
#removed conversion factor - its just shuffle number

#v5.6.1
#implements sum which greater fractional operon counting

#version 5.7
#reorganizes the way the shuffling, etc is done.
#used annotate_group_5.3
#1 get all of the annot into a list
#2 get a list of truncated fasta files
#3 get all of the annot2 and annot3 into a list
#4 run and score one iteration at a time! 

#v5.8
#change score_seqs sp that I dont need to add fake annotation
#this allows more efficient use of memory
#add max possilbe to hits and then remove it at the end
#uses scan_and_annotate version 

#v6
#uses truncate_fasta_v5 (which does the annot2/3 business up front)
#uses scan and annotate_v2.5 (which expects truncate_fasta_v5)
#uses annotate_group_6 (which can be either run MC or not MC)

#v6.21
#allows single shuffle
#automatically deals with list_of_pwm_10 being NULL
#fixes issue with seq when shuffling

#v6.22
#passes all of the nog names from annotate_group(real hits) to annotate_group
#in order to avoid a bunch of cbind/rbind operations there

#v6.23
#does some COG group stuff!
#minor streamliming
#now tabulates min FDR below each point

#v6.24
#adds option to skip some of the fdr calculation/tabulation far down the list of real hits
#deals with the genes being returned in annotate_group when doing "real_hits"

#v6.25
#added unique call when loading annot to deal with Past Horia's BS

#v6.26 
#handles passing the correction factors from truncated fasta into scan and annotate

#v6.27
#handles returning promoter to gene distances in annotate_group

#v6.28
#calculates background frequencies and returns log-odds transformed data

#v6.3
#removed option for tabulate results to (not enough perf to be worth the hassle)
#added option for log odds scoring (do_log_odds), effectively merging 3.281/3 with 3.282/4
#added option for sequence retrieval - this takes a list of conditions:
#   fdr_range <- retrieve_sequences[[1]]    =c(0.1,0)
#   rank_range <- retrieve_sequences[[2]]   =c(1,1000)
#   dist_range <- retrieve_sequences[[3]]   =c(-1000,1000) 

#v6.31
#fixed error when doing log-odds on just one pwm

#
master_annotate_group_2 <- function(list_of_pwm_35, list_of_pwm_10 = NULL, vector_of_ECF_names, weights_F = NA,
                                    vector_of_genomes_filenames, shuffle_number = 20, shuffle_at_a_time = 1, #inputs for the whole fxn
                                    vector_of_annot_filenames,vector_of_annot2_filenames,vector_of_annot3_filenames, #inputs for the whole fxn
                                    how_much_seq_F = c(-332,1), spacers_F = c(22,23,24), cores_to_use = 1, operon_limit = 50, #inputs for the truncate fasta bit
                                    return_genes = TRUE, do_log_odds = FALSE, retrieve_sequences = NULL){ 

  ###################################################################################################################
  #lets make up some even weights if we are not provided weights
  if (is.na(weights_F[1])){weights_F = rep(1, length=length(list_of_pwm_35))}
  
  ###################################################################################################################  
  #ok, now, lets do the column shuffling
  #write the matrices
  #for pwm35
  column_shuffle_35 <- matrix(NA, shuffle_number,dim(list_of_pwm_35[[1]])[2])
  
  #for pwm10, if it exists. otherwise, create the variable with NULL
  if (!is.null(list_of_pwm_10)){
    column_shuffle_10 <- matrix(NA, shuffle_number, dim(list_of_pwm_10[[1]])[2])
  } else {
    column_shuffle_10 <- NULL
  }
  
  #fill them in
  for (i in 1:shuffle_number){
    
    #for pwm35
    column_shuffle_35[i,] <- sample(1:dim(list_of_pwm_35[[1]])[2],dim(list_of_pwm_35[[1]])[2])
    
    #for pwm10, if it exists
    if (!is.null(list_of_pwm_10)){
      column_shuffle_10[i,] <- sample(1:dim(list_of_pwm_10[[1]])[2],dim(list_of_pwm_10[[1]])[2])
    }
    
  }
  
  ###################################################################################################################
  #1 get annot into a list and get all of the operons
  ###################################################################################################################
  
  print("reading in the files!")
  
  list_of_annot <- lapply(1:length(list_of_pwm_35), function(i) unique(as.data.frame(readr::read_tsv(vector_of_annot_filenames[i], 
                                                                                              col_names = FALSE,col_type = readr::cols()))))
  list_of_operons <- lapply(1:length(list_of_pwm_35), function(i) get_all_operons(list_of_annot[[i]], operon_limit))
  
  ##################################################################################################################
  #2 get a list of truncated fasta files
  ###################################################################################################################
  
  list_of_truncated_fasta <- parallel::mclapply(1:length(list_of_pwm_35),
                                                function(i) get_truncated_fasta_file(vector_of_ECF_names[i], vector_of_genomes_filenames[i], 
                                                                                     list_of_operons[[i]],
                                                                                     list_of_annot[[i]],vector_of_annot2_filenames[i],
                                                                                     vector_of_annot3_filenames[i],  how_much_seq_F),
                                                mc.cores = cores_to_use, mc.silent = TRUE)

  ##################################################################################################################
  #2a transform matrices to log odds if applicable
  ###################################################################################################################
  if (do_log_odds){
  
    list_of_background_distributions <- lapply(list_of_truncated_fasta, get_background_frequencies)
    
    list_of_pwm_35 <- correct_pwms_for_background_frequencies(list_of_pwm_35, list_of_background_distributions)
    
    if (!is.null(list_of_pwm_10)){
      list_of_pwm_10 <- correct_pwms_for_background_frequencies(list_of_pwm_10, list_of_background_distributions)
    }
    
  }
  ###################################################################################################################
  #3 do the real hits
  ###################################################################################################################
  #ok, now we need to generate one real annot_hits_summary
  print("doing the real_hits")
  
  #do the real hits
  real_hits <- annotate_group_2(list_of_pwm_35, list_of_pwm_10, list_of_truncated_fasta, list_of_operons, spacers_F,
                                 shuffled_order_35 = NA, shuffled_order_10 = NA, cores_to_use = cores_to_use)
  
  #get these scores up front so we can do the ordering
  original_scores <- score_matrix(real_hits[[1]], weights_F)
  original_scores_order <- order(original_scores, decreasing = TRUE)
  
  #get the first three cols to feed into annotate group
  F3andA <- real_hits[[1]][original_scores_order,1:3]
  
  #and now reorder/subset
  real_hits[[1]] <- real_hits[[1]][original_scores_order,] #this is the scores and annotations 
  real_hits[[2]] <- real_hits[[2]][original_scores_order] #this is is.in.operon, vector of numbers from 0 to 1
  real_hits[[3]] <- real_hits[[3]][original_scores_order] #median.dist.to.ecf
  real_hits[[4]] <- real_hits[[4]][original_scores_order] #hits_cogs 
  real_hits[[5]] <- real_hits[[5]][original_scores_order,] #hit_gene
  real_hits[[6]] <- real_hits[[6]][original_scores_order,] #hit_location 
  real_hits[[7]] <- real_hits[[7]][original_scores_order,] #hit_spacer
  original_scores <- original_scores[original_scores_order]
  
  ###################################################################################################################
  #4 ok, now to shuffle
  ###################################################################################################################
  print("doing the shuffle")
  #the iterative part
  
  #split the shuffle up
  split_shuffle <- vector("list", shuffle_number/shuffle_at_a_time)
  for (i in 1:length(split_shuffle)){split_shuffle[[i]] <- seq((1+(shuffle_at_a_time*(i-1))),(shuffle_at_a_time*i))}
  
  list_of_number_better <- parallel::mclapply(1:length(split_shuffle), function(j) do_an_iteration(list_of_pwm_35, column_shuffle_35[split_shuffle[[j]],,drop = FALSE],
                                                                                                   list_of_pwm_10, column_shuffle_10[split_shuffle[[j]],,drop = FALSE],
                                                                                                   list_of_truncated_fasta, list_of_operons, spacers_F, 
                                                                                                   weights_F, original_scores, F3andA), 
                                              mc.cores = cores_to_use, mc.silent = TRUE)
  
  #sum all of the thigns in the list then break out the 
  number_better <- Reduce("+", list_of_number_better)

  ###################################################################################################################
  #5 fdr conversion
  ###################################################################################################################
  
  #calculate effective shuffles
  if (!is.null(list_of_pwm_10)){eff_shuffles <- (shuffle_number*shuffle_at_a_time)}
  if (is.null(list_of_pwm_10)){eff_shuffles <- shuffle_number}
  
  #how numch more shuffled we have than actual
  fdrs <- (number_better/eff_shuffles)/(sum_which_greater(original_scores, real_hits[[2]], original_scores)+1)
  
  #convert the FDRs to min FDRs above
  fdrs <- convert_fdrs(fdrs, original_scores)
  
  ###################################################################################################################
  #5a gene sequences, if applicable
  ###################################################################################################################
  
  #now to return the gene sequences...
  genome_sequences <- NULL
  
  if (!is.null(retrieve_sequences)){
    
    genome_sequences <- get_hit_sequences(retrieve_sequences, real_hits, fdrs, list_of_truncated_fasta)
    
  }
  
  ###################################################################################################################
  #6 return statements
  ###################################################################################################################
  
  #convert to scores
  #that we can use to order and return
  if (return_genes){
    return(list(cbind(real_hits[[3]],paste(round(number_better), round(original_scores,2)), fdrs, real_hits[[1]], real_hits[[4]]), real_hits[[5]], real_hits[[6]], genome_sequences))
  } else {
    return(cbind(real_hits[[3]],paste(round(number_better), round(original_scores,2)), fdrs, real_hits[[1]], real_hits[[4]], genome_sequences))
  }
  
}

##########################################################################################################
#convert_fdrs
convert_fdrs <- function(fdrs, original_scores){
  
  #get the order of the original scores
  order_to_go_in <- order(original_scores, decreasing = TRUE)
  
  #now, going in reverse order of those scores
  for (i in (length(fdrs)-1):1){
    if (fdrs[order_to_go_in[i+1]] < fdrs[order_to_go_in[i]]){
      fdrs[order_to_go_in[i]] <- fdrs[order_to_go_in[i+1]]
    }
  }
  
  #
  return(fdrs)
}

#############################################################################################################
#get background frequencies
get_background_frequencies <- function(output_of_truncated_fasta){
  
  #
  bck_counts <- table(unlist(output_of_truncated_fasta[[1]]))

  at_counts <- sum(as.numeric(bck_counts[match(c("a","t"),names(bck_counts))]), na.rm = TRUE)
  cg_counts <- sum(as.numeric(bck_counts[match(c("c","g"),names(bck_counts))]), na.rm = TRUE)
  
  #
  return(log(c(at_counts,cg_counts,cg_counts,at_counts)/(2*sum(at_counts, cg_counts)),2))

}

#correct for background frequencies (aka turn to log odds scoring)
correct_pwms_for_background_frequencies <- function(list_of_pwms, list_of_log_backgrounds){
  
  #
  for (k in 1:length(list_of_pwms)){
    for (j in 1:4){
      list_of_pwms[[k]][j,] <- list_of_pwms[[k]][j,]-list_of_log_backgrounds[[k]][j]
    }
  }
  
  #
  return(list_of_pwms)
  
}



