#!/bin/bash -e

#SBATCH --job-name=cellphy_hapcon_IUPAC
#SBATCH -A uoa03626
#SBATCH --cpus-per-task=10
#SBATCH --mem=4G
#SBATCH --time=01-00:00:00
#SBATCH --output=data/3_cellphy_log/hapcon_IUPAC/%j_hapcon.out

module load R

export single_ccoal_output_folder="$PWD/data/2_ccoal"
export simphy_folder="$PWD/data/1_simphy"

start_time=$(date +%s)

process_hapcon_folder() {
    local iteration_identifier="_pop99731_sites74578r8"
    local iteration_folder="$single_ccoal_output_folder/_pop99731_sites74578r8"
    local hapconsensus_file="$iteration_folder/full_haplotypes_dir/full_hap.0001_consensus.fasta.IUPAC"
    local hapcon_output_prefix="$iteration_folder/hapcon_IUPAC"
    
    cellphy.sh RAXML --msa "$hapconsensus_file" --model GT16+FO+E --seed 2 --threads 10 --force perf_threads --prefix "$hapcon_output_prefix"
}

execute_r_script() {
    local iteration_identifier="_pop99731_sites74578r8"
    local iteration_folder="$single_ccoal_output_folder/_pop99731_sites74578r8"
    local species_tree_file="$simphy_folder/_pop99731_sites74578r8/s_tree.trees"
    local hapcon_besttree_file="$iteration_folder/hapcon_IUPAC.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
    local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

    Rscript scripts/3_cellphy_scripts/hap_consensus_v2/cellphy_hapcon.R "$species_tree_file" "$hapcon_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export hapcon_besttree_file
export population_size
export replicate
export iteration_identifier

log_file="data/3_cellphy_log/hapcon_IUPAC/master_log.txt"
echo "Job started at $(date)" >> "$log_file"

mv "data/3_cellphy_log/hapcon_IUPAC/${SLURM_JOB_ID}_hapcon.out" "data/3_cellphy_log/hapcon_IUPAC/_pop99731_sites74578r8_${SLURM_JOB_ID}_hapcon.out"

process_hapcon_folder "iteration_identifier"

execute_r_script "iteration_identifier"

echo "Job finished at $(date)" >> "$log_file"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
