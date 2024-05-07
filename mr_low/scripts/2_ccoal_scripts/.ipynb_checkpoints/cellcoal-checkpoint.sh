#!/bin/bash -e
#SBATCH --job-name=cellcoal_array
#SBATCH --partition=milan
#SBATCH --array=500-650
#SBATCH --output=logs/2_ccoal_log/%A_%a.out
#SBATCH --error=logs/2_ccoal_log/%A_%a.err
#SBATCH -A uoa03626    
#SBATCH --ntasks 1       
#SBATCH --cpus-per-task 2      
#SBATCH --mem=150G                 
#SBATCH --time=00-01:00:00

module load R
module load Python

iteration_identifier=$(Rscript "scripts/2_ccoal_scripts"/cellcoal.R "$(ls -1 "data/1_simphy" | sed -n ${SLURM_ARRAY_TASK_ID}p)")
iteration_identifier=$(echo "$iteration_identifier" | grep -o '/_pop[0-9]*_sites[0-9]*r[0-9]*' | tr -d '/_')

mv "logs/2_ccoal_log/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" "logs/2_ccoal_log/${SLURM_ARRAY_JOB_ID}_${iteration_identifier}.err"



