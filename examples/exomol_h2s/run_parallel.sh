#!/bin/bash

N_FILES=34
N_TASKS=8
CONFIG_FILE=exomol_h2s.py

mkdir logs
#echo "Running $N_FILES files in parallel with $N_TASKS tasks at a time"

# Loop over all files
for ((i=0; i<N_FILES; i+=N_TASKS)); do
    # Loop over the tasks (as long as there are files left)
    for ((j=0; j<N_TASKS && i+j<N_FILES; j++)); do
        idx_min=$((i+j))
        idx_max=$((i+j+1))
        #echo "Running file-range $idx_min, $idx_max"
        pyROX $CONFIG_FILE -c -pbar --files_range $idx_min $idx_max > logs/range_${idx_min}_${idx_max}.log 2>&1 &
    done
    wait
done
#echo "All tasks completed."