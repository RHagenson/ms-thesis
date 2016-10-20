#!/usr/bin/env Rscript
# This script rruns two functions found in R-defs/ to create plots in parallel
# generate_data_pairs() creates a data.frame from all .prof files found in a recursively
# found in a given profiles/ directory.
# build_plot() takes the rows of this data.frame and uses them as inputs to generate a plot if
# the result is significant

source("R-defs/generate_data_pairs.R")
source("R-defs/build_plot.R")
library("parallel")

now <- format(Sys.Date(), format="%d-%m-%y")

# Read in commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# The variables changed by CLI
now <- args[2]
profilesDir <- paste(args[1], now, sep="/")

number=250000
figsDir=paste("figs", sep="/")  # Does not need 'now' appended to it, because the tree is built off profileDir
outputDir=paste("outputs", sep="/")  # Same as figsDir above
pValCut=0.05

# Remove figs/ and outputs/
# unlink(paste(figsDir, now, sep="/"), recursive = TRUE, force = TRUE)
# unlink(paste(outputDir, now, sep="/"), recursive = TRUE, force = TRUE)

data_pairs <- generate_data_pairs(
  profilesDir = profilesDir,
  number = number,
  figsDir = figsDir,
  outputDir = outputDir,
  pValCut = pValCut
)

mcmapply(
  build_plot,
  filename = data_pairs$filenameVector,
  number = data_pairs$numberVector,
  profileDir = data_pairs$profileDirVector,
  figsDir = data_pairs$figsDirVector,
  outputDir = data_pairs$outputDirVector,
  pValCut = data_pairs$pValCutVector,
  mc.cores = detectCores() - 1
)
