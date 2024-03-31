#!/bin/bash -e

#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00-01:00:00
#SBATCH --output=data/3_cellphy_log/vcf/consensus_logs/%j_vcf.out
#SBATCH --job-name=cellphy_vcf

export single_ccoal_output_folder="$PWD/data/2_ccoal"
export simphy_folder="$PWD/data/1_simphy"

module load R
module load VCFtools
module load Python/3.10.5-gimkl-2022a

start_time=$(date +%s)

modify_zzz() {
  local initial_vcf="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf.0001"
  local modified_vcf="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf_modified"
  sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' $initial_vcf > $modified_vcf
}

remove_outgcell() {
  local modified_vcf="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf_modified"
  local output_vcf="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf_no_outgcell"
  vcftools --vcf $modified_vcf --remove-indv outgcell --recode --recode-INFO-all --out $output_vcf
}


process_vcf_consensus() {
  local vcf_no_outgcell="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf_no_outgcell.recode.vcf"
  cmd_python="python scripts/3_cellphy_scripts/vcf_bulkseq/vcf_consensus/combine_vcf_sample_weighted_PL_G10.py $vcf_no_outgcell"
  eval "$cmd_python"
}

cleanup_files() {
  local vcf_no_outgcell="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf_no_outgcell.recode.vcf"
  local vcf_modified="$single_ccoal_output_folder/_pop99551_sites84817r8/vcf_dir/vcf_modified"
  
  if [ -f "$vcf_no_outgcell" ]; then
    rm -f "$vcf_no_outgcell"
   echo "Deleted $vcf_no_outgcell"
  fi
  
  if [ -f "$vcf_modified" ]; then
    rm -f "$vcf_modified"
    echo "Deleted $vcf_modified"
  fi
}

mv "data/3_cellphy_log/vcf/consensus_logs/${SLURM_JOB_ID}_vcf.out" "data/3_cellphy_log/vcf/consensus_logs/_pop99551_sites84817r8_${SLURM_JOB_ID}_vcf.out"

modify_zzz 
remove_outgcell 
process_vcf_consensus 
cleanup_files

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"



