#!/bin/bash
#SBATCH --job-name=hm_int
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --array=1-5
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-bigmem
#SBATCH --mem=350G
#SBATCH --time=11:59:00

singularity exec \
  /acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif \
  Rscript Step2_Int.R full.list_f_perm_0${SLURM_ARRAY_TASK_ID}


