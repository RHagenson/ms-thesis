#!/usr/bin/env Rscript
# This script rruns two functions found in R-defs/ to create plots in parallel
# generate_data_pairs() creates a data.frame from all .prof files found in a recursively
# found in a given profiles/ directory.
# build_plot() takes the rows of this data.frame and uses them as inputs to generate a plot if
# the result is significant

source("R-defs/correction.R")
library("optparse")

# Define the CLI arguments
option_list = list(
  make_option(c("-d", "--date"), 
              type="character", 
              default=format(Sys.Date(), format="%d-%m-%y"), 
              help="Date in profiles to process [form = DD-MM-YY]", 
              metavar="character"),
  make_option(c("-t", "--cancerType"), 
              type="character", 
              default=NULL, 
              help="The cancer type, e.g. BRCA, to process [default= %default]", 
              metavar="character")
);

# Parse the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Make profiles, data, and cancerType required CLI arguments
if (is.null(opt$date)){
  print_help(opt_parser)
  stop("-d/--date is a required argument.n", call.=FALSE)
}
if (is.null(opt$cancerType)){
  print_help(opt_parser)
  stop("-t/--cancerType is a required argument.n", call.=FALSE)
}

# The variables changed by CLI arguments
date <- as.character(opt$date)
cancerType <- as.character(opt$cancerType)

adjusted_frames <- correction(date=date, cancerType=cancerType)

# Create a directory for the signficant results
outDir <- paste("./outputs", date, "p-adjusted", cancerType, sep = "/")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Output the significant results
# Write the two data.frames to file
write.csv(adjusted_frames$long, file = paste(outDir, "long_log.csv", sep = "/"))
write.csv(adjusted_frames$short, file = paste(outDir, "short_log.csv", sep = "/"))

# Redirect output to file and STOUT
sink(paste(outDir, "adjResults.txt", sep = "/"), split = TRUE)

# Output the significant results for long files
print(paste("The number of significant long values is:", as.character(sum(adjusted_frames$long$pValAdj < pValCutoff))))
print("Value(s):")
print(as.character(adjusted_frames$long$pVal[which(adjusted_frames$long$pValAdj < pValCutoff)]))
print("Related isoform:")
print(as.character(adjusted_frames$long$isoName[which(adjusted_frames$long$pValAdj < pValCutoff)]))
print("Direction of significance:")
print(as.character(adjusted_frames$long$pValDir[which(adjusted_frames$long$pValAdj < pValCutoff)]))

# Output the significant results for short files
print(paste("The number of significant short values is:", as.character(sum(adjusted_frames$short$pValAdj < pValCutoff))))
print("Value(s):")
print(as.character(adjusted_frames$short$pVal[which(adjusted_frames$short$pValAdj < pValCutoff)]))
print("Related isoform:")
print(as.character(adjusted_frames$short$isoName[which(adjusted_frames$short$pValAdj < pValCutoff)]))
print("Direction of significance:")
print(as.character(adjusted_frames$short$pValDir[which(adjusted_frames$short$pValAdj < pValCutoff)]))
