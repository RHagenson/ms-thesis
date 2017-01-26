#!/usr/bin/env Rscript
# This script takes in the Disorder...Results.tsv file and 
# outputs a file with equal row lengths with the following format
# gene [tab] type [space] <mutation> [comma] <mut> [tab] next_type ...
# where [...] denotes a sep character and <...> a value in the table
#
# Written by: Ryan Hagenson
# Email: rhagenson@unomaha.edu

library(optparse)
source("R-defs/extract_mutations.R")

# Define the CLI arguments
option_list = list(
  make_option(c("-d", "--disorder"), 
              type="character", 
              default=NULL, 
              help="Location of the Disorder Results", 
              metavar="character")
)

# Parse the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Defines required CLI arguments
if (is.null(opt$disorder)){
  print_help(opt_parser)
  stop("-d/--disorder is a required argument.n", call.=FALSE)
}

# Store CLI arguments in easier-to-use variables
disorder_table <- read.delim("processing/disorderPositions/DisorderPositionResults.tsv", stringsAsFactors = FALSE)
# disorder_table <- read.delim(opt$disorder, stringsAsFactors = FALSE)

# Define internal variables
new_frame <- data.frame(Gene=unique(disorder_table$Gene), 
                        Types=NA,
                        Isoforms=NA,
                        Mutations=NA)

# Build new_frame$Types and new_frame$Isoform entries with uniqueness
for (gene in new_frame$Gene) {
  # Find what indices that gene is found in original table
  gene_indices <- which(disorder_table$Gene == gene)
  
  
  # Gather cancer types and isoform numbers into data.frame
  temp_frame <- data.frame(Types=disorder_table[gene_indices, "Type"],
                           Numbers=disorder_table[gene_indices, "Isoform"])
  
  # Remove duplicate rows
  temp_frame <- temp_frame[!duplicated(temp_frame), ]
  
  # Set new_frame to unique values
  new_frame$Types[which(new_frame$Gene == gene)] = paste0(temp_frame$Types, collapse = ",")
  new_frame$Isoforms[which(new_frame$Gene == gene)] = paste0(temp_frame$Numbers, collapse = ",")
}

# Loop through rows of new_frame and output formatted results
for (index in 1:length(new_frame$Gene)) {
  gene <- new_frame[index, "Gene"]
  
  # Build cancer type and isoform number vectors
  cancer_vector <- as.list(strsplit(new_frame[index, "Types"], ",")[[1]])
  iso_num_vector <- as.list(strsplit(new_frame[index, "Isoforms"], ",")[[1]])
  
  # Set/Reset mutations_list to NULL
  mutations_list <- c()
  
  # Loop through indexes to retrieve mutations in 'A251G' format
  for (sec_index in 1:length(cancer_vector)) {
    mutations_list[sec_index] <- paste(extract_mutations(cancer_type = cancer_vector[sec_index], 
                                                         gene = gene,
                                                         isoform_num = iso_num_vector[sec_index]),
                                   collapse = ",")
    
  }
  
  # Collapse the multiple mutations_list entries into a single string
  new_frame[index, "Mutations"] <- paste(mutations_list, collapse = ";")
  
  output_lines <- vector(length = length(new_frame$Gene))
  
  # Prepare to output to file
  for (index in 1:length(new_frame$Gene)) {
    # Start output line with Gene
    next_line <- new_frame[index, "Gene"]
    
    # Extract Types and Mutations lists
    Types <- as.list(strsplit(new_frame[index, "Types"], split = ",")[[1]])
    Muts <- as.list(strsplit(new_frame[index, "Mutations"], split = ";")[[1]])
    
    # Next loop in parallel over Types and Mutations continually adding to output
    for (sec_index in 1:length(Types)) {
      # Add next type
      next_line <- paste(next_line, Types[sec_index], sep="\t")
      
      # Add next set of mutations
      next_line <- paste(next_line, Muts[sec_index], sep=" ")
    }
    
    # Add line to output vector
    output_lines[index] <- next_line
  }
  
  # Output to file
  write(output_lines, "/tmp/test_out_boih")
}
