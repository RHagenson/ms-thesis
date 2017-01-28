#!/usr/bin/env Rscript
# This script runs the hypergeometric enrichment analysis 
# as presented in Ramakanth's tutorial
# File positions are put in at the command line
#
# Written by: Ryan Hagenson
# Email: rhagenson@unomaha.edu

source("dario/GOUtilities.R")
library(optparse)

# Define the CLI arguments
option_list = list(
  make_option(c("-g", "--go"), 
              type="character", 
              default=NULL, 
              help="Location of the Gene Ontology file", 
              metavar="character"),
  make_option(c("-a", "--annotation"), 
              type="character", 
              default=NULL, 
              help="Location of the annotation file", 
              metavar="character"),
  make_option(c("-t", "--term"), 
              type="character", 
              default=NULL, 
              help="Desired GO term: BP, MF, or CC", 
              metavar="character"),
  make_option(c("-s", "--subset"), 
              type="character", 
              default=NULL, 
              help="A list of files, delimited only by commas, with each containing a gene identifiers (each on a new line)", 
              metavar="character"),
  make_option(c("-b", "--background"), 
              type="character", 
              default=NULL, 
              help="A file that contains background gene identifiers (each on a new line)", 
              metavar="character"),
  make_option(c("-o", "--output"), 
              type="character", 
              default=NULL, 
              help="Output location to be prefixed to tables (default table name is 'table_<subset-file>_<term>.tsv')", 
              metavar="character")
);

# Parse the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Defines required CLI arguments
if (is.null(opt$go)){
  print_help(opt_parser)
  stop("-g/--go is a required argument.n", call.=FALSE)
}
if (is.null(opt$annotation)){
  print_help(opt_parser)
  stop("-a/--annotation is a required argument.n", call.=FALSE)
}
if (is.null(opt$term)){
  print_help(opt_parser)
  stop("-t/--terms is a required argument.n", call.=FALSE)
}
if (is.null(opt$subset)){
  print_help(opt_parser)
  stop("-s/--subset is a required argument.n", call.=FALSE)
}
if (is.null(opt$background)){
  print_help(opt_parser)
  stop("-b/--background is a required argument.n", call.=FALSE)
}

# Set default output to current directory
if (is.null(opt$output)) {
  opt$output = as.character("./")
}

# Store CLI arguments in easier-to-use variables
go_file <- as.character(opt$go)
ann_file <- as.character(opt$annotation)

subset_files_list <- as.list(strsplit(opt$subset, ",")[[1]])

print(subset_files_list)

out_loc <- if(endsWith(as.character(opt$output), "/")) {
              as.character(opt$output) 
            } else {
              paste0(as.character(opt$output), "/")
            }

cat("Loading in the background profile from: ", opt$background, "\n")
background_frame <- read.table(opt$background, stringsAsFactors = F)

term <- if(opt$term == "BP") {
    "biological_process"
  } else if (opt$term == "MF") {
    "molecular_function" 
  } else if (opt$term == "CC") {
    "cellular_component"
  } else {
      print_help(opt_parser)
      stop("Invalid -t/--term argument.n", call.=FALSE)
  }

# Value not changed by CLI
pValCutoff <- 0.05


##########################
### Begin the pipeline ###
##########################
# Steps 1-3 are done online or via python

# Step 4, Build the Gene Ontology Graph
cat("Creating the Gene Ontology graph from: ", go_file, "\n")
GOGraph <- buildGOGraphs(goFile=go_file, 
                         regulatoryLinks=FALSE)

# Step 5, Build the Gene Ontology Dictionary
cat("Creating the Gene Ontology dictionary from: ", go_file, "\n")
GODict <- createGODict(goFile=go_file)

# Step 6, Arrange the GO terms of the filtered annotation file to 
#   promote further analysis
cat("Creating the annotation list from: ", ann_file, "\n")
AnnList <- createAnnList(annFile=ann_file)

# Out of order because it only needs to be run once and takes a long time
# Step 8, Convert protein centric annotation list to GO term centric list
cat("Creating a protein-centric annotation list.\n")
termCentricAnn <- getTermCentricAnn(annotations=AnnList)

###
### Begin looping through the provided subset files and run the remainder steps
###

for(entry in subset_files_list) {
  cat("Processing: ", entry, "\n")
  # Create the subset list
  subset_frame <- read.table(entry, stringsAsFactors = F)
  output_file <- paste0(out_loc,
                        "table_", 
                        sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", 
                            "\\2", 
                            entry), 
                        "_", opt$term,
                        ".tsv")
  
  # Step 7, Perform Hypergeometric test for Enrichment Analysis
  cat("Performing the hypergeometric test for enrichment analysis.\n")
  enrichment <- enrichmentAnalysis(subset=subset_frame$V1, 
                                   allProteins=background_frame$V1, 
                                   annotations=AnnList,
                                   GOGraph=GOGraph[[term]], 
                                   fdr=TRUE, 
                                   underrepresented=FALSE)

  # Step 9, Create the enrichment table to aid visualization
  cat("Outputting the enrichment table at: ", output_file, "\n")
  createEnrichmentTable(subset=subset_frame$V1, 
                        termCentricAnn=termCentricAnn,
                        enrichment=enrichment,
                        GOGraph=GOGraph[[term]],
                        GODict=GODict,
                        pvalue=pValCutoff,
                        outfile=output_file)
}
