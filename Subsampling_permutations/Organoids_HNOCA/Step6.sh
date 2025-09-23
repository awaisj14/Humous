#!/bin/bash
#SBATCH --job-name=hm_step8_perm_vs_ref
#SBATCH --output=logs/step8_perm_%A_%a.out
#SBATCH --error=logs/step8_perm_%A_%a.err
#SBATCH --array=1-5
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared-cpu
#SBATCH --mem=64G
#SBATCH --time=11:59:00

set -euo pipefail
mkdir -p logs

PERM="${SLURM_ARRAY_TASK_ID}"
REF="lands_OG/L_MR_O.rds"
SIF="/acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif"

# Bases to process (PNG + compare) — run in each tree
BASES=(
  "/home/users/j/javed/humous_reviews/HNOCA/"
)

for BASE in "${BASES[@]}"; do
echo "[INFO] Working in BASE: ${BASE} (perm ${PERM})"
singularity exec "${SIF}" bash -lc "
    set -euo pipefail
    cd '${BASE}'
    echo '[INFO] Generating PNGs for perm ${PERM}...'
    Rscript Step6.R ${PERM}
    echo '[INFO] Comparing perm ${PERM} to reference...'
    REF_LMR_PATH='${REF}' Step8_corr.R ${PERM}
  "
done

echo "[OK] All done for perm ${PERM}"
