#!/bin/bash
#SBATCH --job-name=hm_int
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --array=1-4
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-bigmem
#SBATCH --mem=300G
#SBATCH --time=08:00:00


singularity exec \
  /acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif \
  Rscript Step4.R $SLURM_ARRAY_TASK_ID

