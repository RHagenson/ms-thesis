#!/bin/bash

for name in $(find /home/ryan/Thesis/disorderCancer/data/profiles/ -type f -exec basename {} \;)
do
    echo $name
    Rscript ./build-single-dist-plot.R  $name 100000
done
