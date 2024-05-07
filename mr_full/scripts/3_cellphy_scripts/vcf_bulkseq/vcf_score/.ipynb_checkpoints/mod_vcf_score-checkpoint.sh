#!/bin/bash -e
#SBATCH --cpus-per-task=55
#SBATCH --mem=30G
#SBATCH --time=03-00:00:00
#SBATCH --output=3_cellphy_slurmout/vcf/%j_vcf.out
#SBATCH --job-name=cellphy_vcf

export single_ccoal_output_folder="$PWD/2_ccoal_folder"
export simphy_folder="$PWD/1_simphy_folder"

module load R
module load Python/3.10.5-gimkl-2022a

start_time=$(date +%s)

process_vcf_folder() {
local iteration_identifier="$1"
local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
local no_outgcell_consensus_file="$iteration_folder/vcf_dir/vcf_no_outgcell_consensus.vcf"
local vcf_output_prefix="$iteration_folder/vcf"
local species_tree_file="$simphy_folder/${iteration_identifier}/1/s_tree.trees"
local vcf_besttree_file="$iteration_folder/vcf.raxml.bestTree"
local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

cellphy.sh RAXML --msa "$no_outgcell_consensus_file" --model GT16+FO --seed 2 --threads 55 --force perf_threads --prefix "$vcf_output_prefix"

Rscript cellphy_scripts/vcfconsensus/cellphy_vcf.R "$species_tree_file" "$vcf_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export vcf_besttree_file
export population_size
export replicate

mv "3_cellphy_slurmout/vcf/${SLURM_JOB_ID}_vcf.out" "3_cellphy_slurmout/vcf/${iteration_identifier}_vcf.out"

process_vcf_folder "$iteration_identifier"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"




