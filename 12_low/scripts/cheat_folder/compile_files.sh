#!/bin/bash

source_dir1="data/1_simphy"
source_dir2="data/2_ccoal"
destination_dir="data/compiled"

for iteration_folder in "$source_dir1"/*; do
  if [ -d "$iteration_folder" ]; then
    corresponding_iteration_folder="$source_dir2/$(basename "$iteration_folder")"

    if [ -d "$corresponding_iteration_folder" ]; then
      destination_iteration_folder="$destination_dir/$(basename "$iteration_folder")"
      mkdir -p "$destination_iteration_folder/trees"
      mkdir -p "$destination_iteration_folder/sequences"
      mkdir -p "$destination_iteration_folder/logs"
      mkdir -p "$destination_iteration_folder/parameters"

      cp -f "$iteration_folder/s_tree.trees" "$destination_iteration_folder/trees/s_tree.trees"
      cp -f "$iteration_folder"/g_trees*.trees "$destination_iteration_folder/trees/"

      if [ -f "$corresponding_iteration_folder/hapcon.raxml.bestTree" ]; then
        cp -f "$corresponding_iteration_folder/hapcon.raxml.bestTree" "$destination_iteration_folder/trees/hapcon.raxml.bestTree"
      fi

      if [ -f "$corresponding_iteration_folder/single.raxml.bestTree" ]; then
        cp -f "$corresponding_iteration_folder/single.raxml.bestTree" "$destination_iteration_folder/trees/single.raxml.bestTree"
      fi
      
      if [ -f "$corresponding_iteration_folder/vcf.raxml.bestTree" ]; then
        cp -f "$corresponding_iteration_folder/vcf.raxml.bestTree" "$destination_iteration_folder/trees/$corresponding_iteration_folder/vcf.raxml.bestTree"
      fi
        
      full_hap_file="$corresponding_iteration_folder/full_haplotypes_dir/full_hap.0001"
      if [ -f "$full_hap_file" ]; then
        cp -f "$full_hap_file" "$destination_iteration_folder/sequences/full_hap.0001"
      fi

      consensus_fasta_file="$corresponding_iteration_folder/full_haplotypes_dir/full_hap.0001_consensus.fasta"
      if [ -f "$consensus_fasta_file" ]; then
        cp -f "$consensus_fasta_file" "$destination_iteration_folder/sequences/haplotype_consensus.fasta"
      fi
      
    cp -f "$corresponding_iteration_folder/log" "$destination_iteration_folder/logs/"
      cp -f "$corresponding_iteration_folder/hapcon.raxml.log" "$destination_iteration_folder/logs/"
      cp -f "$corresponding_iteration_folder/single.raxml.log" "$destination_iteration_folder/logs/"
      cp -f "$corresponding_iteration_folder/vcf.raxml.log" "$destination_iteration_folder/logs/"
      
      base_name_no_r="$(basename "$iteration_folder" | sed 's/r[0-9]*$//')"

      simphy_conf_file="$iteration_folder/$base_name_no_r.conf"
      if [ -f "$simphy_conf_file" ]; then
        cp -f "$simphy_conf_file" "$destination_iteration_folder/parameters/$(basename "$iteration_folder").conf"
      fi

      simphy_params_file="$iteration_folder/$base_name_no_r.params"
      if [ -f "$simphy_params_file" ]; then
        cp -f "$simphy_params_file" "$destination_iteration_folder/parameters/$(basename "$iteration_folder").params"
      fi
      
      echo "simphy_conf_file: $simphy_conf_file"
      echo "simphy_params_file: $simphy_params_file"
      
      
    fi
  fi
done

echo "SUCCESS"
