#function that takes in a list of pwm_scans, annotates them, and then does the cog thing...
#cuttoff is how far upstream of the start site we want to look. -1000 is there to keep back compatbiliity
#operon tells us how to handle (if at all) operon structures
  #none means no operon structures are taken into consideration
  #first only means that operons are taken into consideration only if the hit is in the first gene
  #all means the entire operon is taken into consideration regardless of where the gene hit is
  #an operon is defined as:

#version 4.1 changes
  #changed summarize_hits to discard cutoff first, then do the rest of the concatenation
  #added class name
  #self -21 where no score (by definition)
  #changed some code to take advantage of the fact that the new version of summarize hits returns only 1 of each NOG by default
  #preinitialized matrices where possible (unlcear if it helped :/)
  #skip the first thing
  #added option to skip the other 

#version 5 huge changes
#completley rewritten to use scan_pwm3, annotate hits 5.2
#use match and removes a loop from every operoation
#use readr file reads

#version 5.11
#can deal with first ecf not having any hits

#version 5.2
#uses scan and annotate
#requieres truncated fasta v3

#version 5.3
#is able to take as an input a list of truncated fasta
#uses scan_and_annotate 2.2

#version 5.4
#represents the merge point for annotate_group and annotate_group_fast
#is able to take as an input a list operons
#uses scan_and_annotate 2.4 (which has better speed for operons)
#eliminated go_fast
#precomputes all the hits

#v6
#uses scan_and_annotate_v2.5, which depends on truncate_fasta_v5
#many fewer inputs, since theyre all included in other things

#v6.2
#gains the ability to do multiple scan and annotate shuffles at the same time
#as implemented in scan and annotate 2.7
#additional input, which is the shuffled order
#remove annot_hit_summary. we can just use the appropriate call!

#v6.21
#this version calculates and returns 
#1. distance from the ECF
#2. spacers that give optimal score
#also, remove annotation from hits location...
#return of additional_information??

#version6.22
#got rid of some checks, ensure compatibility with data_frame output

#v6.23
#can deal with list_of_pwm_10 == NULL

#v6.24
#can do multiple single pwms at a time

#v6.25
#has functionality to take in set of all nogs, etc.

#v6.251
#adds the ability to return information about the strand of of the ECF and the genes

#v6.252
#removed additional_information flag, now uses first_three_cols p/a as flag 
#will also return COG class for all genes

#v6.253
#will now also return as part of real hits (ie when first_three_cols is not present) the names of the genes in the organism
#that make up the results - this will be useful for assessing HGT and also looking up genes for people

#v6.254
#handles returning gene distances for each gene for the real hits

#v6.255
#handles returning spacer lengths as well
#slightly changes the way that median distance to the ECF is calculated


#
annotate_group_2 <- function(list_of_pwm_35, list_of_pwm_10, list_of_truncated_fasta, operon = "yes", spacers, 
                             shuffled_order_35 = NA, shuffled_order_10 = NA, cores_to_use = 1, 
                             first_three_cols_and_annot = NULL){
  
  #############################################################################################################################
  #if the operon is nit a list make it a list
  if (class(operon) != "list"){operon <- as.list(rep(operon, length = length(list_of_pwm_35)))}
  
  #############################################################################################################################
  #we dont have a shuffled order...
  was_na <- is.na(shuffled_order_35[1])
  
  if (was_na){
    
    #for pwm35
    shuffled_order_35 <- matrix(1:dim(list_of_pwm_35[[1]])[2], 1, dim(list_of_pwm_35[[1]])[2])
    
    #for pwm10, if it exists
    if (!is.null(list_of_pwm_10)){
      shuffled_order_10 <- matrix(1:dim(list_of_pwm_10[[1]])[2], 1, dim(list_of_pwm_10[[1]])[2])
    }
  }
  
  #############################################################################################################################
  #precompute the shuffled ECF lists as a list of lists
  
  #for pwm35
  list_of_shuffled_pwm_35 <- list()

  for (i in 1:length(list_of_pwm_35)){

    list_of_shuffled_pwm_35[[i]] <- lapply(1:dim(shuffled_order_35)[1], function(j) list_of_pwm_35[[i]][,shuffled_order_35[j,]])

  }
  
  #for pwm10, if it exists
  if (!is.null(list_of_pwm_10)){
    
    #make the list 
    list_of_shuffled_pwm_10 <- list()
    
    #
    for (i in 1:length(list_of_pwm_10)){
      
      list_of_shuffled_pwm_10[[i]] <- lapply(1:dim(shuffled_order_10)[1], function(j) list_of_pwm_10[[i]][,shuffled_order_10[j,]])
    
    }
    
  } else {
    
    #make a list of NA (we cannot pass a list of NULL)
    list_of_shuffled_pwm_10 <- rep(NA, length(list_of_pwm_35))
    
  }
  
  #############################################################################################################################
  #precompute all of the hits
  list_of_hits <- parallel::mclapply(1:length(list_of_pwm_35), function (i) scan_and_annotate(list_of_shuffled_pwm_35[[i]], 
                                                                                              list_of_shuffled_pwm_10[[i]], 
                                                                                              list_of_truncated_fasta[[i]], 
                                                                                              spacer = spacers, operon[[i]]), 
                                     mc.cores = cores_to_use, mc.silent = TRUE)
  
  #############################################################################################################################
  #now lets annotate each of the things one by one
  to_return <- vector("list", length(list_of_hits[[1]]))
  
  #for the number of replicates
  for (a in 1:length(list_of_hits[[1]])){
    
    #start here 
    start_ecf <-  1
    
    #if we have those first three cols
    if (!is.null(first_three_cols_and_annot)){
      
      #set up the matrices and get on with it
      hits_scores <- cbind(first_three_cols_and_annot, 
                           matrix(NA, dim(first_three_cols_and_annot)[1], length(list_of_pwm_35)))
      
      hit_location <-  matrix(NA, dim(first_three_cols_and_annot)[1], length(list_of_pwm_35))
      
    }
    
    #if we dont have the first_three_cols_and_annot, then we need to do the initialization as per usual
    if (is.null(first_three_cols_and_annot)){
    
      #initial thing
      #preallocate the correct number of columns
      hits_annot <- list_of_hits[[start_ecf]][[a]][,13]
      
      hits_cogs <- list_of_hits[[start_ecf]][[a]][,12]
      
      hits_scores <- cbind(list_of_hits[[start_ecf]][[a]][,c(9,4,14)],matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (start_ecf-1)),
                           list_of_hits[[start_ecf]][[a]][,3], matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (length(list_of_pwm_35[-1])-start_ecf+1)))
      
      hit_location <- cbind(matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (start_ecf-1)),
                            list_of_hits[[start_ecf]][[a]][,7],matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (length(list_of_pwm_35[-1])-start_ecf+1)))
      
      hit_spacer <- cbind(matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (start_ecf-1)),
                            list_of_hits[[start_ecf]][[a]][,4],matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (length(list_of_pwm_35[-1])-start_ecf+1)))
    
      hit_distance <- cbind(matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (start_ecf-1)),
                             list_of_hits[[start_ecf]][[a]][,10],matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (length(list_of_pwm_35[-1])-start_ecf+1)))
      
      hit_gene <- cbind(matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (start_ecf-1)),
                            list_of_hits[[start_ecf]][[a]][,8],matrix(NA, dim(list_of_hits[[start_ecf]][[a]])[1], (length(list_of_pwm_35[-1])-start_ecf+1)))
      
      start_ecf <-  2
      
    }
    
    #for each thing,
    for (i in start_ecf:length(list_of_pwm_35)){
      
      #get the correspondance between matched and unmatched
      matched_cogs <- match(list_of_hits[[i]][[a]][,9], hits_scores[,1])
      
      #for all of the hits that match whats already in the thing
      if (length(which(!is.na(matched_cogs))) > 0){
        
        #here we need both the hits we are working with and where in the hits_scores they
        #should go
        hits_of_int <- which(!is.na(matched_cogs))
        positons_of_int <- matched_cogs[hits_of_int]
        
        #since the matrix column number is preallocated,
        hits_scores[positons_of_int,(3+i)] <- list_of_hits[[i]][[a]][hits_of_int,3]
        hit_location[positons_of_int,i] <- list_of_hits[[i]][[a]][hits_of_int,7]
        
        #we only do this for the real, when we dont have F3andA
        if (is.null(first_three_cols_and_annot)){
        
          hit_spacer[positons_of_int,i] <- list_of_hits[[i]][[a]][hits_of_int,4]
          hit_distance[positons_of_int,i] <- list_of_hits[[i]][[a]][hits_of_int,10]
          hit_gene[positons_of_int,i] <- list_of_hits[[i]][[a]][hits_of_int,8]
          
        }
      }
      
      #for all of the hits that dont match whats already in the thing
      if (length(which(is.na(matched_cogs))) > 0){
        
        #here we only need the positions in annot_hits_summary
        hits_of_int <- which(is.na(matched_cogs))
        
        #rbind an appropriately sized matrix to the bottom of hits summary
        hits_scores <- rbind(hits_scores, setNames(cbind(list_of_hits[[i]][[a]][hits_of_int, c(9,4,14)], 
                                                         matrix(NA, length(hits_of_int), (i-1)), 
                                                         list_of_hits[[i]][[a]][hits_of_int,3], 
                                                         matrix(NA, length(hits_of_int), (length(list_of_pwm_35) - i))), names(hits_scores)))
        
        hit_location <- rbind(hit_location, setNames(cbind(matrix(NA, length(hits_of_int), (i-1)), 
                                                           list_of_hits[[i]][[a]][hits_of_int,7],
                                                           matrix(NA, length(hits_of_int), (length(list_of_pwm_35) - i))), names(hit_location)))
        
        #we only do this for the real, when we dont have F3andA
        if (is.null(first_three_cols_and_annot)){
        
          hit_spacer <- rbind(hit_spacer, setNames(cbind(matrix(NA, length(hits_of_int), (i-1)), 
                                                         list_of_hits[[i]][[a]][hits_of_int,4],
                                                         matrix(NA, length(hits_of_int), (length(list_of_pwm_35) - i))), names(hit_spacer)))
                
          hit_distance <- rbind(hit_distance, setNames(cbind(matrix(NA, length(hits_of_int), (i-1)), 
                                                             list_of_hits[[i]][[a]][hits_of_int,10],
                                                             matrix(NA, length(hits_of_int), (length(list_of_pwm_35) - i))), names(hit_distance)))
          
          hit_gene <- rbind(hit_gene, setNames(cbind(matrix(NA, length(hits_of_int), (i-1)), 
                                                             list_of_hits[[i]][[a]][hits_of_int,8],
                                                             matrix(NA, length(hits_of_int), (length(list_of_pwm_35) - i))), names(hit_gene)))
          
          #and the annotation
          hits_annot <- c(hits_annot, list_of_hits[[i]][[a]][hits_of_int,13])
          hits_cogs <- c(hits_cogs, list_of_hits[[i]][[a]][hits_of_int,12])
          
        }
      }
    }
    
    ###############################
    ###############################
    ###############################
    if (is.null(first_three_cols_and_annot)){
    
      #here we get the spacer mode and percentage that are mode
      which.is.not.na <- lapply(1:dim(hit_spacer)[1], function (i) which(!is.na(hit_spacer[i,])))
      spacer.mode <- unlist(lapply(1:dim(hit_spacer)[1], function (i) Mode(hit_spacer[i,which.is.not.na[[i]]])))
      percentage.mode <- rowMeans(Reduce("cbind", lapply(1:dim(hit_spacer)[2], function (i) hit_spacer[,i] == spacer.mode)), na.rm=TRUE)
      
      #here we get the median distance to the ECF start site. 
      median_minus <- function(a_numeric_vector){
        
        #which are on the same contig
        same_contig <- which(a_numeric_vector != 10000000)
        normal_ones <- intersect(which(!is.na(a_numeric_vector)), same_contig)
        
        
        #
        return(ifelse(length(normal_ones) == 0, 10000000, median(abs(a_numeric_vector[normal_ones]))))
        
      }
      
      median.dist.to.ecf <- unlist(lapply(1:dim(hit_distance)[1], function (i) median_minus(as.numeric(hit_distance[i,]))))
    
      #here we get whether the gene and the ecf are on the same strand
      fraction.different.strand <- unlist(lapply(1:dim(hit_distance)[1], function (i) length(which(hit_distance[i,] < 0))))/
                                   unlist(lapply(1:dim(hit_distance)[1], function (i) length(which(hit_distance[i,] < 9999999))))
      
      median.dist.to.ecf <- median.dist.to.ecf*ifelse(fraction.different.strand >= 0.5, -1, 1)
      
      #put the spacer mode in hits_score
      hits_scores[,2] <- paste(spacer.mode, round(percentage.mode,2))
    
    }
    
    #here we get the percentage of things with operon at every position and return it as part of the list 
    #get the number that are in operons per row
    #apply a rule (here we do fraction which are operons)
    is.in.operon <- 1-rowMeans(do.call("cbind", lapply(1:dim(hit_location)[2], function (i) hit_location[,i] == -1000)), na.rm=TRUE)
    
    if (is.null(first_three_cols_and_annot)){
      to_return[[a]] <- list(cbind(hits_scores,hits_annot),is.in.operon, median.dist.to.ecf, hits_cogs, hit_gene, hit_location, hit_spacer)
    }
    
    #if we are just doing shuffled
    if (!is.null(first_three_cols_and_annot)){
      to_return[[a]] <- list(hits_scores, is.in.operon)
    }
  }
  
  #the return statements!
  if (!was_na){return(to_return)}
  if (was_na){return(to_return[[1]])}
  
}

#mode function just in case
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


