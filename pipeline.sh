#!/bin/bash

set -v
set -x

##
## This script runs a single cancer type through the entire pipeline
## through python and R script written from the MS thesis of 
## Ryan Hagenson, University of Nebraska -- Omaha
##

# Test if at least the date has been given
if test "$1" == "" ; then
    echo $'\a'Please supply date parameter, e.g:
    echo $0 20-10-16
    echo "Should supply date of run in DD-MM-YY format"
    exit
fi

# Define variables
ORIGINAL_DIR=$(pwd)
DATE=$1
CANCER=$2

# Change into python directory to build profiles
cd python/disorder
python create_csv_profile.py --data ../../../disorderCancer/data/ --cancerTypes "$CANCER" --date "$DATE"

# Change to R directory to build LOGs
cd "$ORIGINAL_DIR"
cd ./R/
Rscript build_all_logs.R --profiles ../../disorderCancer/data/profiles/ --date "$DATE" --cancerType "$CANCER"

# Be sure the user is returned to the original directory
cd $ORIGINAL_DIR
