#!/bin/bash -e

#SBATCH --cpus-per-task=10
#SBATCH --partition=milan
#SBATCH --mem=20G
#SBATCH --time=00-03:00:00
#SBATCH --job-name=vcf_concat
#SBATCH --output=logs/7_vcf_concat_log/%j_concat.out

module load VCFtools

main_folder="data/7_vcf_concat"

subdirs=($(find "$main_folder" -maxdepth 1 -mindepth 1 -type d))

for subdir in "${subdirs[@]}"; do
  vcf_files=($(find "$subdir" -type f -name "vcf_no_outgcell.recode_consensus.vcf.NEW_ERr[0-9]*"))
  if [ ${#vcf_files[@]} -eq 0 ]; then
    echo "No VCF files found with the common identifier in directory: $subdir"
    continue
  fi

  output_file="$subdir/concat.vcf"

  vcf-concat "${vcf_files[@]}" > "$output_file"
done
