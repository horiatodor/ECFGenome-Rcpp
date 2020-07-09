#this function does something very much like clustalW for weighting sequences
#based on a tree
#it first makes an unrooted NJ tree
#then does sum(branch lengths) weighted by how many things share a branch!

#version 2.1 
#if everything has weight 0, then all 1s are returned

#load the needed fucntions
get_weights_clustalw <- function(dist_object){
  
  #if there are only two, then the weights will be the same for the two of them
  if (attr(dist_object, "Size") == 2){return(c(0.5,0.5))}
  
  #first lets build the NJ tree
  nj_tree <- ape::nj(dist_object)
  #then midpoint root it
  nj_tree <- phangorn::midpoint(nj_tree)
  
  #weights
  branch_weights <- rep(0, length = attr(dist_object, "Size"))
  node_depths <- ape::node.depth(nj_tree)
  
  #ok, now for every node terminal node
  for (i in 1:length(branch_weights)){
    
    #the first one is easy
    next_edge <- which(nj_tree$edge[,2] == i)
    next_node <- nj_tree$edge[next_edge,1]
    branch_weights[i] <- branch_weights[i] + nj_tree$edge.length[next_edge]
    
    #now, as long as we are not at the terminal node
    while (next_node %in% nj_tree$edge[,2]){
      
      #how many sequences does this node share?
      number_to <- node_depths[next_node]
      
      #next thing
      next_edge <- which(nj_tree$edge[,2] == next_node)
      next_node <- nj_tree$edge[next_edge,1]
      
      #get the score
      branch_weights[i] <- branch_weights[i] + (nj_tree$edge.length[next_edge]/number_to)
      
    }
    
    
  }
  
  #check that the weights are not 0
  if (sum(branch_weights == 0) == length(branch_weights)){
    
    #if so, make them all 1
    return(rep(1,length = length(branch_weights)))
    
  }
  
  #
  return(branch_weights)
  
}

##############
#this is a function that takes in a matrix of distances and some genome names
#and returns weights based on the distances in that file, using get_weights_clustalw
get_weights <- function(of_int_genomes, big_distance_matrix, missing_default = 80){
  
  #now lets make a matrix
  of_int_sim_matrix <- matrix(missing_default, length(of_int_genomes), length(of_int_genomes))
  for (i in 1:length(of_int_genomes)){of_int_sim_matrix[i,i] <-  100}
  
  #now lets get the relevant bits
  relevant_all_of_it <- matrix(big_distance_matrix[intersect(which(big_distance_matrix[,1] %in% of_int_genomes), 
                                                             which(big_distance_matrix[,2] %in% of_int_genomes)),], ncol=3)
  
  if (length(relevant_all_of_it) != 0){
    
    for (i in 1:dim(relevant_all_of_it)[1]){
      
      location1 <- match(relevant_all_of_it[i,1], of_int_genomes)
      location2 <- match(relevant_all_of_it[i,2], of_int_genomes)
      
      of_int_sim_matrix[location1,location2] <- as.numeric(relevant_all_of_it[i,3])
      of_int_sim_matrix[location2,location1] <- as.numeric(relevant_all_of_it[i,3])
      
    }
    
  }
  
  return(get_weights_clustalw(as.dist(100-of_int_sim_matrix)))
  
}

#this is just a little helper function that translates ECF names into the genome name
get_genome_name <- function(ECF_name){
  
  #return thing
  thing_to_return <- rep("", length=length(ECF_name))
  
  #for each thing
  for (i in 1:length(ECF_name)){
    
    #split it
    splitsies <- strsplit(ECF_name[i], "\\.")[[1]]
    thing_to_return[i] <- paste0(splitsies[1:2], collapse=".")
    
  }
  
  #
  return(thing_to_return)
  
}