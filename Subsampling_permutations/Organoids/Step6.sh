#!/bin/bash
#SBATCH --job-name=hm_step8_perm_vs_ref
#SBATCH --output=logs/step8_perm_%A_%a.out
#SBATCH --error=logs/step8_perm_%A_%a.err
#SBATCH --array=1-5
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared-bigmem
#SBATCH --mem=300G
#SBATCH --time=11:59:00

set -euo pipefail
mkdir -p logs

REF="/home/users/j/javed/Humous_reviews/Org/landsO_OG/L_MR_O.rds"
PERM="$SLURM_ARRAY_TASK_ID"

# permutations_extended
BASE="/home/users/j/javed/Humous_reviews/Org/mergedO_permutations" \
USE_AVG=FALSE \
PERM_IDS="$PERM" \
REF_LMR_PATH="$REF" \
singularity exec /acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif \
  Rscript Step6.R

# permutations
BASE="/home/users/j/javed/Humous_reviews/Org/mergedO_permutations_2k" \
USE_AVG=FALSE \
PERM_IDS="$PERM" \
REF_LMR_PATH="$REF" \
singularity exec /acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif \
  Rscript Step6.R
