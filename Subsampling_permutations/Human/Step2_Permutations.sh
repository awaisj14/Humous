#!/bin/bash
#SBATCH --job-name=hm_int
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=08:00:00
GCC/12.2.0  OpenMPI/4.1.4
module load R/4.2.2   # or your cluster’s R
cd /home/users/j/javed/humous_reviews

# run the Rscript wrapper
./Step2_Permutations.R


