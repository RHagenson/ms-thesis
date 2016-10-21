#!/usr/bin/env Rscript
# This script rruns two functions found in R-defs/ to create plots in parallel
# generate_data_pairs() creates a data.frame from all .prof files found in a recursively
# found in a given profiles/ directory.
# build_plot() takes the rows of this data.frame and uses them as inputs to generate a plot if
# the result is significant

source("R-defs/generate_data_pairs.R")
source("R-defs/generate_log.R")
library("parallel")
library("optparse")

# Define the CLI arguments
option_list = list(
  make_option(c("-f", "--figs"), 
              type="character", 
              default=paste("figs", sep="/"), 
              help="Location of the figs/ directory, where PDFs are output [default= %default]", 
              metavar="character"),
  make_option(c("-d", "--date"), 
              type="character", 
              default=format(Sys.Date(), format="%d-%m-%y"), 
              help="Date in profiles to process [form = DD-MM-YY]", 
              metavar="character"),
  make_option(c("-n", "--number"), 
              type="double", 
              default=250000, 
              help="Number of samples to take for each profile [default= %default]", 
              metavar="number"),
  make_option(c("-p", "--profiles"), 
              type="character", 
              default=NULL, 
              help="Location of profiles/ directory", metavar="character"),
  make_option(c("-o", "--outputs"), 
              type="character", 
              default=paste("outputs", sep="/"), 
              help="Location of the outputs/ directory, where LOGs are output [default= %default]", 
              metavar="character"),
  make_option(c("-c", "--pValueCutoff"), 
              type="double", 
              default=0.05, 
              help="The empicial p-Value cutoff [default= %default]", 
              metavar="number")
);

# Parse the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Make profiles and data required CLI arguments
if (is.null(opt$profiles)){
  print_help(opt_parser)
  stop("-p/--profiles is a required argument.n", call.=FALSE)
}
if (is.null(opt$date)){
  print_help(opt_parser)
  stop("-d/--date is a required argument.n", call.=FALSE)
}

# The variables changed by CLI arguments
now <- opt$date
profilesDir <- opt$profiles
number <- opt$number
figsDir <- opt$figs # Does not need 'now' appended to it, because the tree is built off profileDir
outputDir <- opt$outputs # Same as figsDir above
pValCut <- opt$pValueCutoff

# Remove figs/ and outputs/
# unlink(paste(figsDir, now, sep="/"), recursive = TRUE, force = TRUE)
# unlink(paste(outputDir, now, sep="/"), recursive = TRUE, force = TRUE)

data_pairs <- generate_data_pairs(
  profilesDir = paste(profilesDir, now, sep="/"),
  number = number,
  figsDir = figsDir,
  outputDir = outputDir,
  pValCut = pValCut,
  plot = TRUE
)

mcmapply(
  generate_log,
  filename = data_pairs$filenameVector,
  number = data_pairs$numberVector,
  profileDir = data_pairs$profileDirVector,
  figsDir = data_pairs$figsDirVector,
  outputDir = data_pairs$outputDirVector,
  pValCut = data_pairs$pValCutVector,
  plot = data_pairs$plotVector,
  mc.cores = detectCores() - 1
)
