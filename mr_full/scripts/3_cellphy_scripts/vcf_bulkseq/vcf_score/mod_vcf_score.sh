#!/bin/bash -e
#SBATCH --cpus-per-task=10
#SBATCH --partition=milan
#SBATCH --mem=8G
#SBATCH --time=00-06:00:00
#SBATCH --output=logs/3_cellphy_log/vcf/score_logs/%j_vcf.out
#SBATCH --job-name=cellphy_vcf

export single_ccoal_output_folder="$PWD/data/2_ccoal"
export simphy_folder="$PWD/data/1_simphy"

module load R

start_time=$(date +%s)

process_vcf_folder() {
local iteration_identifier="_pop42676_sites99185r4"
local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
local no_outgcell_consensus_file="$iteration_folder/vcf_dir/vcf_no_outgcell.recode_consensus.vcf.NEW_ER"
local vcf_output_prefix="$iteration_folder/vcf"

cellphy.sh RAXML --msa "$no_outgcell_consensus_file" --model GT16+FO --seed 2 --threads 10 --force perf_threads --prefix "$vcf_output_prefix" --force model_lh_impr #--blopt nr_safe --pat-comp off
}

execute_r_script() {
local iteration_identifier="_pop42676_sites99185r4"
local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
local species_tree_file="$simphy_folder/${iteration_identifier}/s_tree.trees"
local vcf_besttree_file="$iteration_folder/vcf.raxml.bestTree"
local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

Rscript scripts/3_cellphy_scripts/vcf_bulkseq/vcf_score/cellphy_vcf.R "$species_tree_file" "$vcf_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export vcf_besttree_file
export population_size
export replicate
export iteration_identifier 


log_file="logs/3_cellphy_log/vcf/master_log.txt"
echo "Job started at $(date)" >> "$log_file"

mv "logs/3_cellphy_log/vcf/score_logs/${SLURM_JOB_ID}_vcf.out" "logs/3_cellphy_log/vcf/score_logs/_pop42676_sites99185r4_${SLURM_JOB_ID}_vcf.out"

process_vcf_folder 
execute_r_script 

echo "Job finished at $(date)" >> "$log_file"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
