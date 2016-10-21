#!/bin/bash

# Test if at least the date has been given
if test "$1" == "" ; then
    echo $'\a'Please supply date parameter, e.g:
    echo $0 20-10-16
    echo "Should supply date of run in DD-MM-YY format"
    exit
fi

# Define variables
ORIGINAL_DIR=${pwd}
DATE=$1

# Change into python directory to build profiles
cd python/disorder
python create_csv_profile.py --data ../../../disorderCancer/data/ --cancerTypes BRCA --date $DATE

# Change to R directory to build LOGs
cd ../../R/
Rscript build_all_logs.R --profiles ../../disorderCancer/data/profiles/ --date $DATE

# Be sure the user is returned to the original directory
cd $ORIGINAL_DIR
