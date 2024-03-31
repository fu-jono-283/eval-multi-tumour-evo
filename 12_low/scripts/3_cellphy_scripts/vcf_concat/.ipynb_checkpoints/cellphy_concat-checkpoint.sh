#!/bin/bash -e

#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=03-00:00:00
#SBATCH --output=logs/3_cellphy_log/vcf_concat/%j_concat.out
#SBATCH --job-name=cellphy_vcf_concat

export concat_iterations="$PWD/data/7_vcf_concat"
export simphy_folder="$PWD/data/1_simphy"

module load R

start_time=$(date +%s)

process_concat_folder() {
    local iteration_identifier="$iteration_identifier"
    local iteration_folder="$concat_iterations/$iteration_identifier"
    local concat_file="$iteration_folder/concat.vcf"
    local concat_output_prefix="$iteration_folder/vcf_concat"
       
    cellphy.sh RAXML --msa "$concat_file" --model GT16+FO --seed 2 --threads $cpus_per_task --force perf_threads --prefix "$concat_output_prefix" --force model_lh_impr
}

execute_r_script() {
    local iteration_identifier="$iteration_identifier"
    local iteration_folder="$concat_iterations/$iteration_identifier"
    local species_tree_file="$iteration_folder/s_tree.trees"
    local concat_besttree_file="$iteration_folder/vcf_concat.raxml.bestTree"
    local population_size=$(echo "$iteration_folder" | grep -oE '_pop([0-9]+)_sites[0-9]' | grep -oE '[0-9]+' | head -1)
    
    Rscript scripts/3_cellphy_scripts/vcf_concat/cellphy_concat.R "$species_tree_file" "$concat_besttree_file" "$population_size" 
}

export species_tree_file
export concat_besttree_file
export population_size


mv "logs/3_cellphy_log/vcf_concat/${SLURM_JOB_ID}_concat.out" "logs/3_cellphy_log/vcf_concat/$iteration_identifier_${SLURM_JOB_ID}_concat.out"

process_concat_folder "$iteration_identifier"
execute_r_script "$iteration_identifier"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
