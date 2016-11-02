#!/bin/bash

###
### This script takes one input: date in form DD-MM-YY and determines which
### cancer types have not been processed and runs through through
### pipeline.sh, which processes everything for a single cancer
###

# Error handling
set -e

# User-input, DATE should be in form DD-MM-YY
DATE=$1
PRO_PATH="./R/outputs/$DATE/"

# Gather all the cancer type identifiers
CANCERS=$(find /home/rhagenson/disorderCancer/data/allMuts/ -mindepth 1 -maxdepth 1 -type f | grep -oP '(\w+)(?=_)')
PROCESSED=$(find $PRO_PATH -mindepth 1 -maxdepth 1 -type d | grep -oP '(\w+)$' || echo "NONE")

echo "Cancer types that have already processed: "
echo $PROCESSED

# Get index of value
#echo ${CANCERS[@]/BRCA//} | cut -d/ -f1 | wc -w |tr -d ' '

echo "Looping through processed to ensure no duplicate runs"
TEMP=("${CANCERS[@]}")
for entry in $PROCESSED; do
    TEMP=$(echo ${TEMP[@]/$entry/})
done

# Inform user of what is happening
echo "Entries in cancer types:" 
echo "$TEMP"

# Run the pipeline for each cancer type in turn and append STDOUT and STDERR
for cancer in ${TEMP[@]}; do
    echo "Now processing $cancer"
    source ./pipeline.sh $DATE $cancer 2>> "$DATE".err >> "$DATE".out &
    echo "PID is $!"
    wait $!
done
