#!/bin/bash -e

#This script basically runs Cellphy for the sc-SEQ data... and then it updates the iteration's CSV file with the tree score. It is executed by the Rscript 'pre...'

#SBATCH --job-name=cellphy_single
#SBATCH --partition=large
#SBATCH -A uoa03626 
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=03-00:00:00
#SBATCH --output=logs/3_cellphy_log/single/%j_single.out

module load R

export single_ccoal_output_folder="$PWD/data/2_ccoal"
export simphy_folder="$PWD/data/1_simphy"

start_time=$(date +%s)

process_single_folder() {
    local iteration_identifier="$iteration_identifier"
    local iteration_folder="$single_ccoal_output_folder/$iteration_identifier"
    local haplotype_file="$iteration_folder/full_haplotypes_dir/full_hap.0001"
    local single_output_prefix="$iteration_folder/single" 
    
    cellphy.sh RAXML --msa "$haplotype_file" --model GT16+FO+E --seed 2 --threads $cpus_per_task --force perf_threads --prefix "$single_output_prefix" 
}

execute_r_script() {
    local iteration_identifier="$iteration_identifier"
    local iteration_folder="$single_ccoal_output_folder/$iteration_identifier"
    local gene_tree_file=$(find "$simphy_folder/$iteration_identifier" -type f -name 'g_trees*.trees' -print -quit)
    local single_besttree_file="$iteration_folder/single.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')  

Rscript scripts/3_cellphy_scripts/single_cseq/cellphy_single.R "$gene_tree_file" "$single_besttree_file" "$population_size" "$replicate" 
}

export gene_tree_file
export single_besttree_file
export population_size
export replicate
export iteration_identifier

log_file="logs/3_cellphy_log/single/master_log.txt"
echo "Job started at $(date)" >> "$log_file"

mv "logs/3_cellphy_log/single/${SLURM_JOB_ID}_single.out" "logs/3_cellphy_log/single/$iteration_identifier_${SLURM_JOB_ID}_single.out"

process_single_folder "iteration_identifier"

execute_r_script "iteration_identifier"

echo "Job finished at $(date)" >> "$log_file"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"