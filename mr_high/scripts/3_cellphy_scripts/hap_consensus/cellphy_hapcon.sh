#!/bin/bash -e

#SBATCH --job-name=cellphy_hapcon
#SBATCH -A uoa03626
#SBATCH --partition=milan
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=01-00:00:00
#SBATCH --output=logs/3_cellphy_log/hapcon/%j_hapcon.out

module load R

export single_ccoal_output_folder="$PWD/data/2_ccoal"
export simphy_folder="$PWD/data/1_simphy"

start_time=$(date +%s)

process_hapcon_folder() {
    local iteration_identifier="$iteration_identifier"
    local iteration_folder="$single_ccoal_output_folder/$iteration_identifier"
    local hapconsensus_file="$iteration_folder/full_haplotypes_dir/full_hap.0001_consensus.fasta"
    local hapcon_output_prefix="$iteration_folder/hapcon"
    
    cellphy.sh RAXML --msa "$hapconsensus_file" --model GT16+FO+E --seed 2 --threads $cpus_per_task --force perf_threads --prefix "$hapcon_output_prefix" --redo
}

execute_r_script() {
    local iteration_identifier="$iteration_identifier"
    local iteration_folder="$single_ccoal_output_folder/$iteration_identifier"
    local species_tree_file="$simphy_folder/$iteration_identifier/s_tree.trees"
    local hapcon_besttree_file="$iteration_folder/hapcon.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

    Rscript scripts/3_cellphy_scripts/hap_consensus/cellphy_hapcon.R "$species_tree_file" "$hapcon_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export hapcon_besttree_file
export population_size
export replicate
export iteration_identifier

log_file="logs/3_cellphy_log/hapcon/master_log.txt"
echo "Job started at $(date)" >> "$log_file"

mv "logs/3_cellphy_log/hapcon/${SLURM_JOB_ID}_hapcon.out" "logs/3_cellphy_log/hapcon/$iteration_identifier_${SLURM_JOB_ID}_hapcon.out"

process_hapcon_folder "iteration_identifier"

execute_r_script "iteration_identifier"

echo "Job finished at $(date)" >> "$log_file"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
