#!/bin/bash
#SBATCH --job-name=hm_int
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --array=1-5
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-cpu
#SBATCH --mem=128G
#SBATCH --time=08:00:00


# --- Singularity container setup ---
SIMG="/acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif"       # <-- path to your Singularity R image

# Avoid thread oversubscription inside container
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# --- Config paths ---
R_SCRIPT="Step4.R"   # your adapted Step4 script
LANDS_DIR="lands"               # folder containing IntH_L_perm_%02d.rds
OUT_OBJ_DIR="LandsH"             # where LandsS_H_* will be written

mkdir -p logs "${OUT_OBJ_DIR}"

# --- Discover input files and pick the one for this array task ---
mapfile -t FILES < <(ls -1 ${LANDS_DIR}/IntH_L_perm_*.rds 2>/dev/null | sort)
N_FILES=${#FILES[@]}
if [[ ${N_FILES} -eq 0 ]]; then
  echo "No input files found in ${LANDS_DIR}/IntH_L_perm_*.rds" >&2
  exit 1
fi

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ ${IDX} -lt 0 || ${IDX} -ge ${N_FILES} ]]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} out of range (have ${N_FILES} files)" >&2
  exit 1
fi

FILE="${FILES[$IDX]}"

# Extract permutation ID from filename (_perm_XX.rds)
if [[ "$FILE" =~ _perm_([0-9]+)\.rds$ ]]; then
  PADDED="${BASH_REMATCH[1]}"
else
  echo "Could not extract perm_id from filename: $FILE" >&2
  exit 1
fi
# Strip leading zeros to get decimal int
PERM_ID=$((10#$PADDED))

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: FILE=${FILE} -> perm_id=${PERM_ID}"
echo "Running inside Singularity: ${SIMG}"

# --- Run Step4 script inside container ---
singularity exec \
  --bind $PWD:/mnt \
  "${SIMG}" \
  Rscript /mnt/${R_SCRIPT} "${PERM_ID}" /mnt/${LANDS_DIR} /mnt/${OUT_OBJ_DIR}
