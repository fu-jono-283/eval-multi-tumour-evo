#!/bin/bash -e

#SBATCH -A uoa03626
#SBATCH --job-name=cutoff
#SBATCH --output=scripts/cheat_folder/ncount/logs/output.log  
#SBATCH --error=scripts/cheat_folder/ncount/logs/errors.log      
#SBATCH --ntasks 1       
#SBATCH --cpus-per-task 2      
#SBATCH --mem=4G          
#SBATCH --time=00-01:30:00

module load Python/3.10.5-gimkl-2022a

run_consensus_cutoff() {
  python scripts/cheat_folder/ncount/cutoff/cutoff_value.py "$1" "$2" "$3"
}

export single_ccoal_output_folder="$PWD/data/2_ccoal"
cutoff_values=(0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40)

full_hap_file="$single_ccoal_output_folder/$iteration_identifier/full_haplotypes_dir/full_hap.0001"

if [ -e "$full_hap_file" ]; then
  for cutoff in "${cutoff_values[@]}"; do
    output_folder=$(dirname "$full_hap_file")
    output_filename="ncount_cutoff_$cutoff.fasta"
    output_path="$output_folder/$output_filename"
    run_consensus_cutoff "$cutoff" "$full_hap_file" "$output_path" "$iteration_identifier"
  done
fi

echo "Full Hap File: $full_hap_file"
echo "Output Folder: $output_folder"
echo "Output Filename: $output_filename"
echo "Output Path: $output_path"
echo "Processing iteration $iteration_identifier"


