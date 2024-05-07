#!/bin/bash -e
#SBATCH --job-name=consensus
#SBATCH --partition=large
#SBATCH --output=data/2_ccoal_slurmout/%j_con.out
#SBATCH --error=data/2_ccoal_slurmout/%j_con.err
#SBATCH -A uoa03626    
#SBATCH --ntasks 1       
#SBATCH --cpus-per-task 2      
#SBATCH --mem=4G                 
#SBATCH --export NONE
#SBATCH --time=00-05:00:00

module load R
module load Python

exported_iteration_identifier="$iteration_identifier"

output_parent_folder="data/2_ccoal"

iteration_folder="$output_parent_folder/$exported_iteration_identifier"

update_consensus_sequence() {
    full_hap_file="$1"
    cmd_python="python scripts/2_ccoal_scripts/consensus_cutoff.py \"$full_hap_file\""
    echo "Running: $cmd_python"
    eval "$cmd_python"
}

echo "Processing iteration folder: $iteration_folder"

full_hap_file="$iteration_folder/full_haplotypes_dir/full_hap.0001"

if [ -f "$full_hap_file" ]; then
    echo "Updating consensus sequence for: $full_hap_file"
    update_consensus_sequence "$full_hap_file" 
    echo "Consensus sequence updated."
else
    echo "Full hap file not found for iteration: $iteration_folder"
fi
