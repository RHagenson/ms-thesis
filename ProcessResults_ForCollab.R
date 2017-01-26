#!/usr/bin/env Rscript
# This script takes in the Disorder...Results.tsv file and 
# outputs a file with equal row lengths with the following format
# gene [tab] type [space] <mutation> [comma] <mut> [tab] next_type ...
# where [...] denotes a sep character and <...> a value in the table
#
# Written by: Ryan Hagenson
# Email: rhagenson@unomaha.edu

library(optparse)

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
new_frame <- data.frame(Gene=unique(disorder_table$Gene), Types=length(unique(disorder_table$Gene)))

for (gene in new_frame$Gene) {
  # Find what indices that gene is found in original table
  gene_indices <- which(disorder_table$Gene == gene)
  
  # Get cancer types
  cancers <- disorder_table[gene_indices, "Type"]
  
  # Set new_frame at index to cancers
  new_frame$Types[which(new_frame$Gene == gene)] = paste(unique(cancers), collapse=',')
}
