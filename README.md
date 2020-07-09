# ECFGenome-Rcpp

Author: Horia Todor [horia.todor@gmail.com]

This code performs phylogenetic footprinting on the basis of a motif consisting of 1 or 2 PWMs separated by a variable spacer. 
It requieres genome and annotation files from proGenomes (Mende et al., 2017). Example files are provided for 10 pseudomonas genomes.

## How to use this code

The main function is called `master_annotate_group_2`. It takes the following inputs:

`list_of_pwm_35` is a list of PWMs

`list_of_pwm_10 = NULL`, is an optional list of PWMs for the second motif

`vector_of_ECF_names` is a vector containing the gene names of ECFs (or other regulator) being queried. It is used to determine auto-regulation. 

`weights_F = NA` is a vector of weights that represents the phylogenetic relationship between the genomes queried

`vector_of_genomes_filenames` is a vector containing the filenames of the fasta files containing the genomes of interest.

`shuffle_number = 20` number of time to shuffle the PWMs to generate the random distribution.

`shuffle_at_a_time = 1` number of shuffles done at a time. For bipartite motifs, all shuffles of one motif are evaluated with all shuffles of the other motif. So 
`shuffle_number = 20, shuffle_at_a_time = 5` evaluates 5x5 4 times for a total of 100 randomized motifs. 

`vector_of_annot_filenames` is a vector containing the filenames of the annot files (see example data)

`vector_of_annot2_filenames` is a vector containing the filenames of the annot2 files (see example data)

`vector_of_annot3_filenames` is a vector containing the filenames of the annot3 files (see example data)

`how_much_seq_F = c(-332,1) designates the amount of sequence to examine upstream of each annotated ORF 

`spacers_F = c(22,23,24)` if using a bipartite motif, this is a vector of spacers between the leading edge of each motif

`cores_to_use = 1` is the number of CPU cores to use. Only for linux/unix/macOS systems due to the use of a 'clone' call to do memory efficient multi-cpu utilization

`operon_limit = 50` not implemented. 

`return_genes = TRUE` option to return genes

`do_log_odds = FALSE` is the option to score PWMs as log-odds. Due to normalization in other parts of the pipeline, this is redundant. 

`retrieve_sequences = NULL` is an option to retrieve sequences associated with the results

