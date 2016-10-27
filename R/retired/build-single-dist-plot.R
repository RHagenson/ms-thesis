#!/usr/bin/env Rscript
# This script takes a disorder mutation profile created by parsing TCGA data
# and generates a .pdf with the normal disorder plot marked by the observed disorder score

args <- commandArgs(trailingOnly = TRUE)
#print(args)

# The parameters that will eventual change
filename <- args[1] # "MUC16.001.long"  # The profile being processed
N = as.numeric(args[2]) # 100000  # The number of samples to take

profileDir = args[3]
figsDir = args[4]  # The location of where figures should be written, do not ending forget '/'
pValCut = as.numeric(args[5])

source("R-defs/build-single-dist-plot.R")

build.plot(filename, N, profileDir, figsDir, pValCut)
