#!/bin/bash -e
#SBATCH --job-name=cellphy_single
#SBATCH --partition=milan
#SBATCH --output=3_cellphy_slurmout/single/%j_single.out
#SBATCH -A uoa03626 
#SBATCH --cpus-per-task=50
#SBATCH --mem=20G
#SBATCH --time=06-00:00:00

module load R

export single_ccoal_output_folder="$PWD/2_ccoal_folder"
export simphy_folder="$PWD/1_simphy_folder"

start_time=$(date +%s)

process_single_folder() {
    local iteration_identifier="$1"
    local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
    local haplotype_file="$iteration_folder/full_haplotypes_dir/full_hap.0001"
    local single_output_prefix="$iteration_folder/single" 
    
    cellphy.sh RAXML --msa "$haplotype_file" --model GT16+FO+E --seed 2 --threads 50 --force perf_threads --prefix "$single_output_prefix"
}

execute_r_script() {
    local iteration_identifier="$1"
    local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
    local gene_tree_file="$simphy_folder/${iteration_identifier}/1/g_trees1.trees"
    local single_besttree_file="$iteration_folder/single.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')  

Rscript cellphy_scripts/single_cseq/cellphy_single.R "$gene_tree_file" "$single_besttree_file" "$population_size" "$replicate" 
}

export gene_tree_file
export single_besttree_file
export population_size
export replicate
export iteration_identifier

log_file="3_cellphy_slurmout/single/master_log.txt"
echo "Job started at $(date)" >> "$log_file"

mv "3_cellphy_slurmout/single/${SLURM_JOB_ID}_single.out" "3_cellphy_slurmout/single/${iteration_identifier}_${SLURM_JOB_ID}_single.out"

process_single_folder "$iteration_identifier"

execute_r_script "$iteration_identifier"

echo "Job finished at $(date)" >> "$log_file"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
