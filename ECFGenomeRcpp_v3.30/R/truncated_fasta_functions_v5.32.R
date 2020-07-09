#this is a function that reads in a fasta file and an annotation file 
#and returns a truncated fasta file and a translation from truncated to full length fasta file

#v 2.1 
#removed some hard coded constants

#v3 returns a list of upstream sequences
#added name protection for ("\t") shit

#v4 removes the third part of the list, where the matched old sequence is stored
#this verison should be used with scan_and_annotate_v2.2.R
#which will not try and access that thing. 

#v5 now does all of the annotation up front here!
#annot2 and annot3

#v5.1 includes the genome location in the name of the list of sequences
#moved the removal of non "+" and "-" annot[,7] outside of the iterating loop
#better annotation throughout the function

#v5.2 makes some minor speed boosts

#v5.3 now does counting and concatenation by NOGS as well, appends to the list. needs the operons as an input

#v5.31 now returns a slightly different output for a slightly different NOG correction

#v5.32 fixes the 0 vs 1 indexing issue

get_truncated_fasta_file <- function(ECF_name, genome_fasta_filenames, operon,
                                     annot, annotation_file2, annotation_file3, how_much_seq = c(-332,1)){
  
  #first lets read in the things that we need
  if (class(annot) == "character"){
  
      annot <- as.data.frame(readr::read_tsv(annot, col_names = FALSE,col_type = readr::cols()))

    } 
  
  #get rid of things in annot that do not make sense!
  to_keep <- which(annot[,7] %in% c("+", "-"))
  annot <- annot[to_keep, ]
  
  #now lets increment annot cols 5 and 6 to make sense in light of the 0 vs 1 indexing
  annot[,5] <- annot[,5] + 1
  annot[,6] <- annot[,6] + 1
  
  #read in the genome
  files_to_score <- seqinr::read.fasta(genome_fasta_filenames)
  
  #initialize the outputs
  truncated_fasta <- NULL
  genome_coordinates <- NULL
  reordered_annot <- NULL
  
  #ok, so we have to go on a by contig basis for speed
  #the goal will be to have translated fasta look just like a fasta file
  #for each contig, we will get
  for (a in 1:length(files_to_score)){
    
    #get the name and length
    contig_name <- names(files_to_score)[a]
    length_of_contig <- length(files_to_score[[a]])
    
    #if the name has a "\t" character in it, remove it
    if (length(grep("\\t", contig_name)) != 0){
      contig_name <- strsplit(contig_name, "\\t")[[1]][1]
    }
    
    #get the relevant_annotation
    hits_of_int <- which(annot[,2] == contig_name)
    
    if (length(hits_of_int) != 0){
   
      #get just the relevant annotation
      relevant_annot <- annot[hits_of_int,]
      
      #change this to an ifelse statement,which, plus and minus
      start_spots <- ifelse(relevant_annot[,7] == "+", 
                            (relevant_annot[,5]+how_much_seq[1]), 
                            (relevant_annot[,6]-how_much_seq[2]))
      
      stop_spots <- ifelse(relevant_annot[,7] == "+", 
                           (relevant_annot[,5]+how_much_seq[2]),
                           (relevant_annot[,6]-how_much_seq[1]))
      
      #now put these unadultatated things into the genome coordinates
      genome_coordinates[[a]] <- cbind(start_spots, stop_spots)
      reordered_annot[[a]] <- relevant_annot
      
      #now to get the number of start Ns
      start_ns <- ifelse(start_spots < 1, 1-start_spots, 0)
      start_spots <- ifelse(start_spots < 1, 1, start_spots)
      
      stop_ns <- ifelse(stop_spots > length_of_contig, stop_spots-length_of_contig, 0)
      stop_spots <- ifelse(stop_spots > length_of_contig, length_of_contig, stop_spots)
      
      #write the sequence
      truncated_fasta <- c(truncated_fasta, lapply(1:length(hits_of_int), function (i) c(rep("n", start_ns[i]), 
                                                                                         files_to_score[[a]][start_spots[i]:stop_spots[i]],
                                                                                         rep("n", stop_ns[i]))))

    }
  }
  
  #rbind relevant things
  reordered_annot <- do.call("rbind", reordered_annot)
  genome_coordinates <- do.call("rbind", genome_coordinates)

  #now match annot2
  #read in the annotation2 file
  if (class(annotation_file2) == "character"){
    
    annot2 <- as.matrix(readr::read_tsv(annotation_file2, col_names = FALSE,col_type = readr::cols()))
    additional_info2 <- annot2[match(reordered_annot[,1], annot2[,2]),c(3,4,11,12,13)]
    
  } else {
    
    additional_info2 <- annotation_file2[match(reordered_annot[,1], annotation_file2[,2]),c(3,4,11,12,13)]
    
  }
  
  #now match annot3
  #get annotation 3
  if (class(annotation_file3) == "character"){
    
    annot3 <- as.matrix(readr::read_tsv(annotation_file3, col_names = FALSE,col_type = readr::cols()))
    additional_info3 <- annot3[match(reordered_annot[,1], annot3[,2]),6]
    
  } else {
    
    additional_info3 <- annotation_file3[match(reordered_annot[,1], annotation_file3[,2]),6]
    
  }
  
  #now do self
  #change to self when its the ECF
  ecf_spot <- which(reordered_annot[,1] == ECF_name)
  if (length(ecf_spot != 0)){additional_info2[ecf_spot,1] <- "SELF"}
  
  #do the counting
  #table of_nogs
  nog_counts <- table(additional_info2[,1])
  
  #set up the thing we will return
  nog_list_and_counts <- rep(0, length(nog_counts))
  
  #do the count by gene
  number_per_gene <- rep(NA, length = dim(operon[[1]])[1])
  
  for (i in 1:length(operon[[2]])){
    
    number_per_gene[operon[[2]][[i]]] <- seq(1:length(operon[[2]][[i]]))
    
  }
  
  #add to number_per_gene all of the other genes
  not_in_yet <- reordered_annot[which(!(reordered_annot[,1] %in% operon[[1]][,1])),1]
  operon[[1]] <- rbind(operon[[1]], cbind(not_in_yet, rep(1,length(not_in_yet))))
  number_per_gene <- c(number_per_gene, rep(1,length(not_in_yet)))
  
  #reorder to match the big thing
  reorder_operons <- match(reordered_annot[,1], operon[[1]][,1])
  number_per_gene <- number_per_gene[reorder_operons]
  
  #match the genes and NOGs
  matched_genes <- match(additional_info2[,1], names(nog_counts))
  
  #add the counts by operon and by NOG
  for (i in 1:length(matched_genes)){
    
    #get all of the genes in the NOG
    nog_list_and_counts[matched_genes[i]] <- nog_list_and_counts[matched_genes[i]]+sum(number_per_gene[i])
  
  }
  
  
  #
  return(list(truncated_fasta, cbind(reordered_annot, additional_info2, additional_info3), 
              how_much_seq, as.data.frame(genome_coordinates), nog_counts, nog_list_and_counts))
  
}

##########################################################################################################