#!/usr/bin/env Rscript
# This script rruns two functions found in R-defs/ to create plots in parallel
# generate_data_pairs() creates a data.frame from all .prof files found in a recursively
# found in a given profiles/ directory.
# build_plot() takes the rows of this data.frame and uses them as inputs to generate a plot if
# the result is significant

source("R-defs/generate_data_pairs.R")
source("R-defs/build_plot.R")
library("parallel")

data_pairs <- generate_data_pairs(profilesDir = "~/Thesis/data/profiles/")

mcmapply(build_plot, 
         filename=data_pairs$filenameVector,
         number=data_pairs$numberVector,
         profileDir=data_pairs$profileDirVector,
         figsDir=data_pairs$figsDirVector,
         pValCut=data_pairs$pValCutVector)