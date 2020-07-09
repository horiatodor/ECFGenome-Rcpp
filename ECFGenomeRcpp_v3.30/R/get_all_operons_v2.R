#this is a function that takes in an annotation file from progenomes and returns a list with all operons > 1, 
#in upstream to downstream order

get_all_operons <- function(annot, distance = 50){
  
  #results function
  results <- list()
  
  #for each contig
  unique_contigs <- unique(annot[,2])
  
  for (ct in 1:length(unique_contigs)){
    
    #get the "relevant annotation" of that contig
    relevant_annot_contig <- which(annot[,2] == unique_contigs[ct])
    
    #in each direction
    #fwd
    fwd_of_int <- intersect(which(annot[,7] == "+"),relevant_annot_contig)
    
    if (length(fwd_of_int) > 1){
      
      #this is the relevant annot for all fwd on this contig
      relevant_annot <- annot[fwd_of_int,]
      
      #sort it by largest to smallest start site
      relevant_annot <- relevant_annot[order(relevant_annot[,5], decreasing = FALSE),]
      
      #gene dist is the distance between the finish of a gene and the start of the gene right behind it
      gene_dist <- c(relevant_annot[-dim(relevant_annot)[1],6]-relevant_annot[-1,5],1000)
      
      #
      genes_of_int <- intersect(which(gene_dist >= -50), which(gene_dist <= 50 ))
      
      #
      i=1
      while (i <= length(genes_of_int)){
        
        current_operon <- relevant_annot[genes_of_int[i],1]  
        
        if ((genes_of_int[i]+1) < length(gene_dist)){
          while (gene_dist[genes_of_int[i]+1] >= -50 & 
                 gene_dist[genes_of_int[i]+1] <= 50){
            
            current_operon <- c(current_operon, relevant_annot[genes_of_int[i]+1,1])
            if ((genes_of_int[i]+1) > length(gene_dist)){break}
            if ((i+1) > length(genes_of_int)){break
              }else {i=i+1}
            
          }
        }
        current_operon <- c(current_operon, relevant_annot[(genes_of_int[i]+1),1] )
        
        results[[length(results)+1]] <- current_operon
        if (i+1 > length(genes_of_int)){break
        }else {i=i+1}
      }
      
    }
    
    #rev
    rev_of_int <- intersect(which(annot[,7] == "-"),relevant_annot_contig)
    
    if (length(rev_of_int) > 1){
      
      #this is the relevant annot for all fwd on this contig
      relevant_annot <- annot[rev_of_int,]
      
      #sort it by largest to smallest start site
      relevant_annot <- relevant_annot[order(relevant_annot[,5], decreasing = TRUE),]
      
      #gene dist is the distance between the finish of a gene and the start of the gene right behind it
      gene_dist <- c(relevant_annot[-dim(relevant_annot)[1],5]-relevant_annot[-1,6],1000)
      
      #
      genes_of_int <- intersect(which(gene_dist >= -50), which(gene_dist <= 50 ))
      
      #
      i=1
      while (i <= length(genes_of_int)){
        
        current_operon <- relevant_annot[genes_of_int[i],1]  
        
        if ((genes_of_int[i]+1) < length(gene_dist)){
          while (gene_dist[genes_of_int[i]+1] >= -50 & 
                 gene_dist[genes_of_int[i]+1] <= 50){
            
            current_operon <- c(current_operon, relevant_annot[genes_of_int[i]+1,1])
            if (genes_of_int[i]+1 > length(gene_dist)){break}
            if (i+1 > length(genes_of_int)){break
            }else {i=i+1}
            
          }
        }
        
        current_operon <- c(current_operon, relevant_annot[(genes_of_int[i]+1),1] )
        
        results[[length(results)+1]] <- current_operon
        if (i+1 > length(genes_of_int)){break
        }else {i=i+1}
      }
      
    }
     
  }
  
  #get unique number for each operon to go with unlist operons
  operon_ids <- NULL
  operon_list <- list()
  start_here <- 1
  
  for (i in 1:length(results)){
    
    operon_ids <- c(operon_ids, rep(i, length(results[[i]])))
    operon_list[[i]] <- seq(start_here,(start_here+length(results[[i]])-1),by=1)
    start_here <- start_here+length(results[[i]])    
  
  }
  
  #return the results
  #return(results)
  return(list(cbind(unlist(results), operon_ids), operon_list))
  
}