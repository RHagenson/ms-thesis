#!/usr/bin/env Rscript
# This script rruns two functions found in R-defs/ to create plots in parallel
# generate_data_pairs() creates a data.frame from all .prof files found in a recursively
# found in a given profiles/ directory.
# build_plot() takes the rows of this data.frame and uses them as inputs to generate a plot if
# the result is significant

source("R-defs/generate_data_pairs.R")
source("R-defs/build_plot.R")
library("parallel")

# Read in commandline arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  print("Please provide the commandline arguments")
  print("Order of arguments: profilesDir, number, figsDir, pValCut")
  print("At least profilesDir is needed")
  print("Example: ~/Thesis/disorderCancer/data/profiles 100000 figs/ 0.05")
  q()
}

len <- as.character(length(args))
# Set defaults of args, at least profilesDir is needed
switch(len,
       "0"={
         print("Please provide the commandline arguments")
         print("Order of arguments: profilesDir, number, figsDir, pValCut")
         print("At least profilesDir is needed")
         print("Example: ~/Thesis/disorderCancer/data/profiles 100000 figs/ 0.05")
         q()
       },
       "1"={
         profilesDir = args[1]
         number=1000000
         figsDir="figs/"
         pValCut=0.05
         },
       "2"={
         profilesDir = args[1]
         number = args[2]
         figsDir="figs/"
         pValCut=0.05
         },
       "3"={
         profilesDir = args[1]
         number = args[2]
         figsDir = args[3]
         pValCut=0.05
         },
       "4"={
         profilesDir = args[1]
         number = args[2]
         figsDir = args[3]
         pValCut=args[4]
         }
       )


data_pairs <- generate_data_pairs(
  profilesDir = profilesDir,
  number = number,
  figsDir = figsDir,
  pValCut = pValCut
)

mcmapply(
  build_plot,
  filename = data_pairs$filenameVector,
  number = data_pairs$numberVector,
  profileDir = data_pairs$profileDirVector,
  figsDir = data_pairs$figsDirVector,
  pValCut = data_pairs$pValCutVector
)
