#!/bin/bash -e
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=01-00:00:00
#SBATCH --output=3_cellphy_slurmout/vcf/%j_vcf.out
#SBATCH --job-name=cellphy_vcf

export single_ccoal_output_folder="$PWD/2_ccoal_folder"
export simphy_folder="$PWD/1_simphy_folder"

module load R
module load VCFtools
module load Python/3.10.5-gimkl-2022a

start_time=$(date +%s)

modify_zzz() {
  local initial_vcf="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf.0001"
  local modified_vcf="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf_modified"
  sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' $initial_vcf > $modified_vcf
}

remove_outgcell() {
  local modified_vcf="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf_modified"
  local output_vcf="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf_no_outgcell"
  vcftools --vcf $modified_vcf --remove-indv $output_vcf --recode --recode-INFO-all --out $output_vcf
}


process_vcf_consensus() {
  local vcf_no_outgcell="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf_no_outgcell.recode.vcf"
  cmd_python="python cellphy_scripts/vcf_bulkseq/vcf_consensus/consensus_vcf.py $vcf_no_outgcell"
  eval "$cmd_python"
}

cleanup_files() {
  local vcf_no_outgcell="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf_no_outgcell.recode.vcf"
  local vcf_modified="$single_ccoal_output_folder/_pop32147_sites46549r7/vcf_dir/vcf_modified"
  
#  if [ -f "$vcf_no_outgcell" ]; then
#    rm -f "$vcf_no_outgcell"
#    echo "Deleted $vcf_no_outgcell"
#  fi
  
#  if [ -f "$vcf_modified" ]; then
#    rm -f "$vcf_modified"
#    echo "Deleted $vcf_modified"
#  fi
}

mv "3_cellphy_slurmout/vcf/${SLURM_JOB_ID}_vcf.out" "3_cellphy_slurmout/vcf/$iteration_identifier_${SLURM_JOB_ID}_vcf.out"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"



