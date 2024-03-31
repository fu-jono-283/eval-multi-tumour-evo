#!/bin/bash -e

#SBATCH -A uoa03626
#SBATCH --job-name=cutoff_nobackup
#SBATCH --output=scripts/cheat_folder/allcount/logs/output.log  
#SBATCH --error=scripts/cheat_folder/allcount/logs/errors.log      
#SBATCH --ntasks 1       
#SBATCH --cpus-per-task 2      
#SBATCH --mem=4G          
#SBATCH --time=00-00:30:00

module load R
module load Python

Rscript scripts/cheat_folder/allcount/allcount/allcount.R
