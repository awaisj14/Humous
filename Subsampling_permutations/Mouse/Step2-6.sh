#!/bin/bash
#SBATCH --job-name=perm_runs
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --array=1-5              # Run 5 permutations
#SBATCH --cpus-per-task=8        # Use 8 cores for integration/training
#SBATCH --partition=shared-cpu
#SBATCH --mem=128G               # Heavy memory usage for Seurat + integration
#SBATCH --time=11:59:00

set -euo pipefail
mkdir -p logs

# Avoid thread oversubscription inside BLAS/OpenMP
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# ----------------------------
# Paths and constants
# ----------------------------
SIF="/acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif"
SCRIPT="Step2-6.R"   # <-- NEW SCRIPT
LIB_MISC="lib_misc.R"
LIB_LANDS="lib_lands.R"
MERGED_RDS="merged_all.rds"
REF_LMR="landsM_OG/L_MR_M_orth.rds"
OUT_ROOT="permutations_all"
BASE_SEED=101

# ----------------------------
# Get current permutation index from SLURM
# ----------------------------
PERM_INDEX="${SLURM_ARRAY_TASK_ID}"

# ----------------------------
# Run the R script for this permutation
# ----------------------------
singularity exec "${SIF}" Rscript "${SCRIPT}" \
  --lib "${LIB_MISC}" \
  --landlib "${LIB_LANDS}" \
  --out_root "${OUT_ROOT}" \
  --base_seed "${BASE_SEED}" \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --merged "${MERGED_RDS}" \
  --ref_lmr "${REF_LMR}" \
  --n_perms 5 \
  --perm "${PERM_INDEX}"
