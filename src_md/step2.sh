#!/bin/bash

for run in 8DFN_full 8DFN_half 7SI9_full 7SI9_half
do
    qsub "./"$run".pbs"
    echo ">>> '$run.pbs' submitted"
done
