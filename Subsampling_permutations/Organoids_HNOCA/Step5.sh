#!/bin/bash
#SBATCH --job-name=hm_int
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --array=1-4
#SBATCH --cpus-per-task=16
#SBATCH --partition=shared-cpu
#SBATCH --mem=64G
#SBATCH --time=08:00:00

module load GCC/9.3.0
module load singularity-compose/0.0.21-Python-3.8.2

cd /home/users/j/javed/humous_reviews/HNOCA

singularity exec \
  /acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif \
  Rscript Step5.R $SLURM_ARRAY_TASK_ID

