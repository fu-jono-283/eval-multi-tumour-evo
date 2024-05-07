#!/bin/bash -e
#SBATCH --cpus-per-task=25
#SBATCH --mem=20G
#SBATCH --time=03-00:00:00
#SBATCH --output=3_cellphy_slurmout/hapcon/%j_hapcon.out
#SBATCH --job-name=cellphy_hapcon

export single_ccoal_output_folder="$PWD/2_ccoal_folder"
export simphy_folder="$PWD/1_simphy_folder"

module load R

start_time=$(date +%s)

process_hapcon_folder() {
    local iteration_identifier="$1"
    local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
    local hapconsensus_file="$iteration_folder/full_haplotypes_dir/full_hap.0001_consensus.fasta"
    local hapcon_output_prefix="$iteration_folder/hapcon"
    
    cellphy.sh RAXML --msa "$hapconsensus_file" --model GT16+FO+E --seed 2 --threads 25 --force perf_threads --prefix "$hapcon_output_prefix"
}

execute_r_script() {
    local iteration_identifier="$1"
    local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
    local species_tree_file="$simphy_folder/${iteration_identifier}/1/s_tree.trees"
    local hapcon_besttree_file="$iteration_folder/hapcon.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

    Rscript cellphy_scripts/hap_consensus/cellphy_hapcon.R "$species_tree_file" "$hapcon_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export hapcon_besttree_file
export population_size
export replicate
export iteration_identifier

log_file="3_cellphy_slurmout/hapcon/master_log.txt"
echo "Job started at $(date)" >> "$log_file"

mv "3_cellphy_slurmout/hapcon/${SLURM_JOB_ID}_hapcon.out" "3_cellphy_slurmout/hapcon/${iteration_identifier}_${SLURM_JOB_ID}_hapcon.out"

process_hapcon_folder "$iteration_identifier"

execute_r_script "$iteration_identifier"

echo "Job finished at $(date)" >> "$log_file"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
