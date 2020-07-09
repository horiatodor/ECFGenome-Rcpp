#ok, so here we are mergin the scan_pwm function and the annotate hits function
#the goal is to get a value for every gene and to go faster
#this will **hopefully ** generate a distribution to assign better values to outlying points

#version 2.2
#uses truncated_fasta_v4 
#dumped the go terms from annot3
#cleaned up the use of the truncated output

#v2.3
#uses slightly faster scoreseq_plus function that integrates a the pwm scanning
#cleaned up some of the spacer stuff

#version2.4
#speed up operon stuff
#got rid of go fast tag

#version 2.5
#excpects input from truncated fasta v5, which already has done the annot2 and annot3 and SELF matching
#much faster operon conversion

#version 2.6
#uses score_seq_plus_v3, which scores pwm35 and pwm10 together 

#version 2.7 is able to take in multiple pwm sequences

#version 2.71 uses scoreseq4.2 to do combinatoric joining of -35 and -10 using scoreseq 4.2
#additionally this version now does genome location 

#version 2.72 does distances more accurately by considering the strand of the gene

#version 2.8 minor speedups, changes the spacer to now be start of motif to start of motif!!!

#version 2.9 uses new version of do operons and summarize that are faster 
#code was refactored a fair bit

#version 2.91 fixes an issue with the way operons were conveyed starting in 2.9

#version 2.92 attempts to speed up summarize function by preemptively identifying the things that dont need to be messed with

#version 2.93 implements efficient single motif searching.
#pwm 10 may be NA

#version 2.94 allows multiple motifs of single PWM to be scanned together

#version 2.95 scans all of the fwd and rev genes together using a newer version 
#of score_seq_plus and score_seq_single

#version 2.96 does not do distance from ECF promoter in cases where there are multiple ecfs 

#version 2.97 uses a slightly more efficient method for getting all of the operons. 

#version 2.971 fixes the distance

#version 2.972 minor speedups

#version 2.98 now expectes input from truncated fasta 5.3, which includes counts per nog
#additionally/crucially, now implenents correction for operon and NOG preponderance when calculating z-scores

#version 2.99
#uses newer version of do_operon for speed!

#version 2.991
#corrects edge effects of multiple COG correction
#uses summarize 2.8

#version 2.992
#now slightly more accurate with returning the locations of hiits
#changed return location to start of motif


#the funtion itself
scan_and_annotate <- function(pwm35, pwm10, output_of_truncated_fasta, spacer = c(22,23,24), operon = NULL){

  #is the pwm35 a list?
  is.pwm.list <- class(pwm35) == "list" 
  
  #is pwm10 NA
  is.pwm10.na <- is.na(pwm10[1])
  
  #how many things?
  number_of_ecfs <- 1

  #deal with the pwms
  if (is.pwm.list){
    
    #if there is a pwm10, then square the number of things
    #otherwise, its just the number of things
    if (!is.pwm10.na){number_of_ecfs <- length(pwm35)^2}
    if (is.pwm10.na){number_of_ecfs <- length(pwm35)}
    
    #pwm35
    rev_pwm35 <- do.call("rbind", lapply(pwm35, rev_comp_pwm))
    pwm35 <- do.call("rbind", pwm35)
    
    #pwm10
    if (!is.pwm10.na){
      rev_pwm10 <- do.call("rbind", lapply(pwm10, rev_comp_pwm))
      pwm10 <- do.call("rbind", pwm10)   
    }
    
  } else {
    
    #jsut rev comp the pwm35
    rev_pwm35 <- rev_comp_pwm(pwm35)
    
    #and the pwm10, if it exists
    if (!is.pwm10.na){
      rev_pwm10 <- rev_comp_pwm(pwm10)
    } 
  }
  
  ######################################################
  #######################################################
  #######################################################
  #ok here we need to do the hits equivalent
  ##############
  ##############
  how_much_seq <- output_of_truncated_fasta[[3]]
  
  #first, lets  
  #scan in the forward direction
  #if the contig is big enough to contain the largest spacer
  #but we only need to scan the things that face forward!
  fwd_genes <- which(output_of_truncated_fasta[[2]][,7] == "+")
  all_scores <- matrix(NA, dim(output_of_truncated_fasta[[2]])[1], number_of_ecfs)
  all_spacers <- matrix(NA, dim(output_of_truncated_fasta[[2]])[1], number_of_ecfs)
  all_dist_from_tss <- matrix(NA, dim(output_of_truncated_fasta[[2]])[1], number_of_ecfs)

  #here we call the normal version of scoreseq_plus_combinatoric 
  if (!is.pwm10.na){

    #now, for each spacer length, add the matrices and then subset them 
    temp_scores <- scoreseq_plus_combinatoric2(pwm35,pwm10,output_of_truncated_fasta[[1]][fwd_genes], spacer =  spacer)
    
    #fill in the relevant spots
    all_scores[fwd_genes,] <- temp_scores[[1]]
    all_spacers[fwd_genes,] <- temp_scores[[3]]
    all_dist_from_tss[fwd_genes,] <- how_much_seq[1]+temp_scores[[2]]-1
  
  }
  
  #here we call the single version of scoreseq_plus_combinatoric 
  if (is.pwm10.na){
    
    #now, for each spacer length, add the matrices and then subset them 
    temp_scores <- scoreseq_plus_single2(pwm35,output_of_truncated_fasta[[1]][fwd_genes])
    
    #fill in the relevant spots
    all_scores[fwd_genes,] <- temp_scores[[1]] 
    all_spacers[fwd_genes,] <- 0
    all_dist_from_tss[fwd_genes,] <-(how_much_seq[1]+temp_scores[[2]]-1)
  
  }
  
  #scan in the reverse direction
  rev_genes <- which(output_of_truncated_fasta[[2]][,7] == "-")

  #
  #here we call the normal version of scoreseq_plus_combinatoric 
  if (!is.pwm10.na){
    
    #now, for each spacer length, add the matrices and then subset them 
    temp_scores <- scoreseq_plus_combinatoric2(rev_pwm10,rev_pwm35, output_of_truncated_fasta[[1]][rev_genes], spacer = spacer)
    
    #fill in the relevant spots
    all_scores[rev_genes,] <- temp_scores[[1]] 
    all_spacers[rev_genes,] <- temp_scores[[3]]
    all_dist_from_tss[rev_genes,] <- how_much_seq[2] + 2 - temp_scores[[2]] - (dim(pwm35)[2]+temp_scores[[3]])

  }
  
  #here we call the single version of scoreseq_plus_combinatoric 
  if (is.pwm10.na){
    
    #now, for each spacer length, add the matrices and then subset them 
    temp_scores <- scoreseq_plus_single2(rev_pwm35, output_of_truncated_fasta[[1]][rev_genes])
    
    #fill in the relevant spots
    all_scores[rev_genes,] <- temp_scores[[1]] 
    all_spacers[rev_genes,] <- 0
    #this is the length of the sequence (plus one) ((how_much_seq[2] - how_much_seq[1] + 2)
    #minus the "front" of the motif (temp_scores[[2]] + (dim(pwm35)[2] - 1)))
    #this gives us the equivalent location in the fwd direction
    #so now its  
    #all_dist_from_tss[rev_genes,] <- how_much_seq[1]+((how_much_seq[2] - how_much_seq[1] + 2) - (temp_scores[[2]] + (dim(pwm35)[2] - 1)))-1
    #and we simplify and get...
    all_dist_from_tss[rev_genes,] <- how_much_seq[2] + 2 - temp_scores[[2]] - dim(pwm35)[2]
    
  }
  
  #
  
  ################################################################################################
  #common thing we will output  
  hits_common <- cbind(rep("fwd",dim(output_of_truncated_fasta[[2]])[1]), #this reminds us that these are all forward sites
                       output_of_truncated_fasta[[2]][,2],#this puts the contig name on each line 
                       output_of_truncated_fasta[[2]][,1], #this is the gene
                       as.matrix(output_of_truncated_fasta[[2]][,13:18])) #this is all that other crap
  
  hits_common[rev_genes,1] <- "rev" 
  hits_common[,5] <- 10000000
  hits_common <- as.data.frame(hits_common, stringsAsFactors = FALSE)
  
  #if necessary, prepare the operon!
  if (!is.null(operon)){
    
    if (class(operon) != "list"){
      
      #get all the operons
      operon <- get_all_operons(output_of_truncated_fasta[[2]][,1:12])
      
    } 

  }
  
  #precompute a few things
  #start sites
  gene_starts <- rep(0, dim(output_of_truncated_fasta[[2]])[1])
  gene_starts[fwd_genes] <- output_of_truncated_fasta[[4]][fwd_genes,1]
  gene_starts[rev_genes] <- output_of_truncated_fasta[[4]][rev_genes,2]
  
  #operon matches
  if (!is.null(operon)){
    matches <- match(operon[[1]][,1],output_of_truncated_fasta[[2]][,1])
    operon_two <- lapply(operon[[2]], function(i) i-1)
  }
  
  #gene starts
  all_gene_starts <- matrix(rep(gene_starts, number_of_ecfs), ncol = number_of_ecfs)
  all_gene_starts[fwd_genes,] <- all_gene_starts[fwd_genes,]-all_dist_from_tss[fwd_genes,]
  all_gene_starts[rev_genes,] <- all_gene_starts[rev_genes,]+all_dist_from_tss[rev_genes,]
  
  #we dont do distance if there are multiple motifs,
  #because that is a shuffle iteration, where it doesnt really matter...
  if (number_of_ecfs == 1){
    if (length(which(hits_common[,4] == "SELF")) != 0){
      ecf_row <- which(hits_common[,4] == "SELF") 
    } else {
      ecf_row <- NA
    }
  } else {
    ecf_row <- NA
  }
  
  #the thing we will output!
  return(lapply(1:number_of_ecfs, function (i) summarize_hits_final_4(hits_common, all_gene_starts[,i], do_operon_Rcpp(all_scores[,i], operon_two, matches), 
                                                                      all_scores[,i], all_dist_from_tss[,i], all_spacers[,i],
                                                                      names(output_of_truncated_fasta[[5]]), ecf_row, output_of_truncated_fasta[[6]])))
   
}

##################################
#function to rev_comp a pwm
###################################
rev_comp_pwm <- function(pwm){
  
  #reverse
  rev_pwm <- pwm[,dim(pwm)[2]:1]
  
  #complement by changing the rownames 
  rownames(rev_pwm) <- seqinr::comp(rownames(pwm))
  
  #reorder the rownames
  rev_pwm <- rev_pwm[order(rownames(rev_pwm)),]
  
  return(rev_pwm)
  
}

