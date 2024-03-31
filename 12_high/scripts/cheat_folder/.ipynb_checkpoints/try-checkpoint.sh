#!/bin/bash

source_dir1="data/1_simphy"
source_dir2="data/2_ccoal"
destination_dir="data/compiled"

# Iterate through all replicate directories in source_dir1
for replicate_dir1 in "$source_dir1"/*; do
  if [ -d "$replicate_dir1" ]; then
    # Extract the replicate directory name without the path
    replicate_name=$(basename "$replicate_dir1")
    
    # Extract the pop_sites part of the replicate name
    pop_sites=$(echo "$replicate_name" | sed -E 's/_pop[^_]+_sites[^_]+//')

    # Create the destination subdirectory
    dest_subdir="$destination_dir/$pop_sites"
    mkdir -p "$dest_subdir"

    done