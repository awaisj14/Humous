#!/bin/bash
#SBATCH --job-name=hm_step5_pngs
#SBATCH --output=logs/step5_%A_%a.out
#SBATCH --error=logs/step5_%A_%a.err
#SBATCH --array=1-5                 # <-- set to number of permutation files you have
#SBATCH --cpus-per-task=16
#SBATCH --partition=shared-cpu
#SBATCH --mem=64G
#SBATCH --time=11:59:00

set -euo pipefail

mkdir -p logs

# --- Singularity container (keep your existing image path) ---

SIMG="/acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif"

# --- Paths inside your project ---
R_SCRIPT="Step5.R"   # <- updated Step 5 script name/location
LANDS_DIR="lands"                      # <- where L_MR_H_perm_*.rds live

# --- Parallelism & threading (inside container) ---
export STEP5_CORES="${SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# --- Discover permutation files and map array index -> perm_id ---
mapfile -t FILES < <(ls -1 "${LANDS_DIR}"/L_MR_H_perm_*.rds 2>/dev/null | sort)
N_FILES=${#FILES[@]}
if [[ ${N_FILES} -eq 0 ]]; then
  echo "No input files found: ${LANDS_DIR}/L_MR_H_perm_*.rds" >&2
  exit 1
fi

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ ${IDX} -lt 0 || ${IDX} -ge ${N_FILES} ]]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} out of range (have ${N_FILES} files)" >&2
  exit 1
fi

FILE="${FILES[$IDX]}"

# Extract 0-padded perm suffix and convert to decimal (handles 01, 1, etc.)
if [[ "$FILE" =~ _perm_([0-9]+)\.rds$ ]]; then
  PADDED="${BASH_REMATCH[1]}"
else
  echo "Could not extract perm_id from filename: $FILE" >&2
  exit 1
fi
PERM_ID=$((10#$PADDED))

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}/${N_FILES} -> ${FILE} (perm_id=${PERM_ID})"
echo "Running: Rscript ${R_SCRIPT} ${PERM_ID} ${LANDS_DIR} (inside ${SIMG})"

# --- Execute inside the container ---
# Bind the working directory to /mnt to ensure consistent paths inside the container.
singularity exec \
  --bind "$PWD":/mnt \
  "${SIMG}" \
  Rscript /mnt/"${R_SCRIPT}" "${PERM_ID}" /mnt/"${LANDS_DIR}"
