#!/bin/bash -e
#SBATCH -A uoa03626
#SBATCH --job-name=cutoff
#SBATCH --output=n_count/logs/job_output.log  
#SBATCH --error=n_count/logs/job_errors.log      
#SBATCH --ntasks 1       
#SBATCH --cpus-per-task 2      
#SBATCH --mem=4G          
#SBATCH --time=00-01:00:00

module load Python/3.10.5-gimkl-2022a

run_consensus_cutoff() {
  python n_count/cutoff_value.py "$1" "$2" "$3"
}

export single_ccoal_output_folder="$PWD/2_ccoal_folder"
cutoff_values=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)

full_hap_file="$single_ccoal_output_folder/$iteration_identifier/full_haplotypes_dir/full_hap.0001"

if [ -e "$full_hap_file" ]; then
  for cutoff in "${cutoff_values[@]}"; do
    output_folder=$(dirname "$full_hap_file")
    output_filename="consensus_cutoff_$cutoff.fasta"
    output_path="$output_folder/$output_filename"
    run_consensus_cutoff "$cutoff" "$full_hap_file" "$output_path" "$iteration_identifier"
  done
fi

echo "Full Hap File: $full_hap_file"
echo "Output Folder: $output_folder"
echo "Output Filename: $output_filename"
echo "Output Path: $output_path"
echo "Processing iteration $iteration_identifier"


