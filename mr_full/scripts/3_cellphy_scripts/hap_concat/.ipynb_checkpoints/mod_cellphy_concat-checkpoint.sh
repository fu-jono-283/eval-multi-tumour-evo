#!/bin/bash -e
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=03-00:00:00
#SBATCH --output=3_cellphy_slurmout/hapcon/%j_hapcon.out
#SBATCH --job-name=cellphy_hapcon

export single_ccoal_output_folder="$PWD/astral/hapcon/iterations"
export simphy_folder="$PWD/1_simphy_folder"

module load R

start_time=$(date +%s)

process_hapcon_folder() {
    local iteration_identifier="$1"
    local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
    local hapconsensus_file="$iteration_folder/merged.fasta"
    local hapcon_output_prefix="$iteration_folder/merged"
    local species_tree_file="$iteration_folder/s_tree.trees"
    local hapcon_besttree_file="$iteration_folder/hapcon.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')
    
    cellphy.sh RAXML --msa "$hapconsensus_file" --model GT16+FO+E --seed 2 --threads 10 --force perf_threads --prefix "$hapcon_output_prefix"
}

execute_r_script() {
    local iteration_identifier="$1"
    local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
    local hapconsensus_file="$iteration_folder/merged.fasta"
    local hapcon_output_prefix="$iteration_folder/merged"
    local species_tree_file="$iteration_folder/s_tree.trees"
    local hapcon_besttree_file="$iteration_folder/hapcon.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

    Rscript cellphy_scripts/hap_consensus/cellphy_hapcon.R "$species_tree_file" "$hapcon_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export hapcon_besttree_file
export population_size
export replicate

mv "3_cellphy_slurmout/hapcon/${SLURM_JOB_ID}_hapcon.out" "3_cellphy_slurmout/hapcon/${iteration_identifier}_hapcon.out"

process_hapcon_folder "$iteration_identifier"
execute_r_script "$iteration_identifier"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
