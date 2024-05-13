#!/bin/bash -e

#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=07-00:00:00
#SBATCH --output=logs/3_cellphy_log/vcf/%j_vcf.out
#SBATCH --job-name=vcf_consensus

export single_ccoal_output_folder="$PWD/data/2_ccoal"
export simphy_folder="$PWD/data/1_simphy"

module load R
module load VCFtools
module load Python/3.10.5-gimkl-2022a

start_time=$(date +%s)

modify_zzz() {
  local initial_vcf="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf.0001"
  local modified_vcf="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf_modified"
  sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' $initial_vcf > $modified_vcf
}

remove_outgcell() {
  local modified_vcf="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf_modified"
  local output_vcf="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf_no_outgcell"
  vcftools --vcf $modified_vcf --remove-indv outgcell --recode --recode-INFO-all --out $output_vcf
}


process_vcf_consensus() {
  local vcf_no_outgcell="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf_no_outgcell.recode.vcf"
  cmd_python="python scripts/3_cellphy_scripts/vcf_consensus/combine_vcf_sample_weighted_PL_G10.py $vcf_no_outgcell"
  eval "$cmd_python"
}

cleanup_files() {
  local vcf_no_outgcell="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf_no_outgcell.recode.vcf"
  local vcf_modified="$single_ccoal_output_folder/_pop23013_sites55238r2/vcf_dir/vcf_modified"
  
  if [ -f "$vcf_no_outgcell" ]; then
    rm -f "$vcf_no_outgcell"
   echo "Deleted $vcf_no_outgcell"
  fi
  
  if [ -f "$vcf_modified" ]; then
    rm -f "$vcf_modified"
    echo "Deleted $vcf_modified"
  fi
}

process_vcf_folder() {
local iteration_identifier="_pop23013_sites55238r2"
local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
local no_outgcell_consensus_file="$iteration_folder/vcf_dir/vcf_no_outgcell.recode_consensus.vcf.NEW_ER"
local vcf_output_prefix="$iteration_folder/vcf"

cellphy.sh RAXML --msa "$no_outgcell_consensus_file" --model GT16+FO --seed 2 --threads 10 --force perf_threads --prefix "$vcf_output_prefix" --force model_lh_impr #--blopt nr_safe --pat-comp off
}

execute_r_script() {
local iteration_identifier="_pop23013_sites55238r2"
local iteration_folder="$single_ccoal_output_folder/${iteration_identifier}"
local species_tree_file="$simphy_folder/${iteration_identifier}/s_tree.trees"
local vcf_besttree_file="$iteration_folder/vcf.raxml.bestTree"
local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]+r[0-9]+' | grep -oE '[0-9]+' | head -1)
local replicate=$(echo "$iteration_folder" | grep -oE 'r[0-9]+$')

Rscript scripts/3_cellphy_scripts/vcf_consensus/cellphy_vcf.R "$species_tree_file" "$vcf_besttree_file" "$population_size" "$replicate"
}

export species_tree_file
export vcf_besttree_file
export population_size
export replicate
export iteration_identifier

mv "logs/3_cellphy_log/vcf/${SLURM_JOB_ID}_vcf.out" "logs/3_cellphy_log/vcf/_pop23013_sites55238r2_${SLURM_JOB_ID}_vcf.out"

modify_zzz 
remove_outgcell 
process_vcf_consensus 
cleanup_files
process_vcf_folder
execute_r_script

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"



