#!/usr/bin/env bash
#SBATCH --job-name=cutrun_full
#SBATCH --output=cutrun_full_%j.out
#SBATCH --error=cutrun_full_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --partition=shared-cpu
#SBATCH --time=11:59:00

set -euo pipefail

# ==============================
# CONFIGURATION
# ==============================
THREADS=16
ROOT_DIR="/srv/beegfs/scratch/users/j/javed/Fastqs_050525/allfastqs/CR_Org"
TRIM_ADAPTERS="/srv/beegfs/scratch/users/j/javed/Fastqs_050525/tools/TruSeq3-PE.fa"
YEAST_INDEX="/srv/beegfs/scratch/users/j/javed/Fastqs_050525/tools/R64-1-1/R64-1-1"
BOWTIE2_INDEX="/srv/beegfs/scratch/users/j/javed/Fastqs_050525/tools/GRCh38_noalt_as/GRCh38_noalt_as"

# --- SEACR Bash wrapper
SEACR_SH="/srv/beegfs/scratch/users/j/javed/Fastqs_050525/tools/SEACR_1.3.sh"

EFFECTIVE_GENOME_SIZE=2913022398
BW_BINSIZE=10
USE_SPIKEIN="yes"

# Keep only IgGCAU (control) and Irf1PCAU (experiment)
SAMPLES=(Irf1PCAU K4CAU)
CONTROL="IgGCAU"

# ==============================
# PREP
# ==============================
cd "$ROOT_DIR"
mkdir -p simple_out/{trimmed,alignment/{bedgraph,bigwig},peakCalling/SEACR}
cd simple_out
echo "[INFO] Workdir: $(pwd)  Threads: $THREADS"

# ==============================
# TRIMMING + MAPPING
# ==============================
trim_and_map() {
  local sample="$1"
  echo "=============================="
  echo "[TRIM+MAP] $sample"
  echo "=============================="

  local raw_R1="../${sample}_R1.fastq.gz"
  local raw_R2="../${sample}_R2.fastq.gz"
  local R1="trimmed/${sample}_R1_trimmed.fastq.gz"
  local R2="trimmed/${sample}_R2_trimmed.fastq.gz"
  local bam="${sample}.bam"
  local filt_bam="${sample}.filtered.bam"
  local spike_bam="${sample}.spikein.bam"

  # ---- TRIMMING ----
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[TRIM] $sample"
    trimmomatic PE -threads "$THREADS" \
      "$raw_R1" "$raw_R2" \
      "$R1" trimmed/${sample}_R1_unpaired.fastq.gz \
      "$R2" trimmed/${sample}_R2_unpaired.fastq.gz \
      ILLUMINACLIP:"$TRIM_ADAPTERS":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
  else
    echo "[SKIP] Trimmed FASTQs exist: $R1, $R2"
  fi

  # ---- HOST ALIGN ----
  if [[ ! -f "$filt_bam" ]]; then
    echo "[ALIGN] $sample -> $bam"
    bowtie2 --dovetail --threads "$THREADS" -x "$BOWTIE2_INDEX" -1 "$R1" -2 "$R2" | \
      sambamba view -S -f bam /dev/stdin > "${bam}.unsorted"
    sambamba sort -t "$THREADS" -o "$bam" "${bam}.unsorted"
    rm -f "${bam}.unsorted"
    sambamba index -t "$THREADS" "$bam"

    echo "[FILTER] remove duplicates; keep mapped/proper pairs"
    sambamba view -f bam -F "proper_pair and not duplicate" -t "$THREADS" "$bam" > "$filt_bam"
    sambamba index -t "$THREADS" "$filt_bam"
  else
    echo "[SKIP] filtered BAM exists: $filt_bam"
  fi

  # ---- SPIKE-IN ALIGN ----
  if [[ ! -f "$spike_bam" ]]; then
    echo "[SPIKEIN] $sample -> yeast"
    bowtie2 --dovetail --threads "$THREADS" -x "$YEAST_INDEX" -1 "$R1" -2 "$R2" | \
      sambamba view -S -f bam /dev/stdin > "${sample}.spikein.unsorted.bam"
    sambamba sort -t "$THREADS" -o "$spike_bam" "${sample}.spikein.unsorted.bam"
    rm -f "${sample}.spikein.unsorted.bam"
    sambamba index -t "$THREADS" "$spike_bam"
  else
    echo "[SKIP] Spike-in BAM exists: $spike_bam"
  fi
}

# Run control + experiment
trim_and_map "$CONTROL"
for s in "${SAMPLES[@]}"; do
  trim_and_map "$s"
done

# ==============================
# SPIKE SCALE + BIGWIGS
# ==============================
get_spike_scale() {
  local sample="$1"
  local s_bam="${sample}.spikein.bam"
  local scale_file="alignment/bedgraph/${sample}.spike.scale.txt"
  if [[ -s "$scale_file" ]]; then cat "$scale_file"; return 0; fi

  local mapped_reads=0
  if [[ -s "$s_bam" ]]; then
    mapped_reads=$(sambamba view -c -F "not unmapped" "$s_bam" || echo 0)
  fi
  local pairs=$(( mapped_reads / 2 ))
  local scale="1"
  if (( pairs > 0 )); then
    scale=$(awk -v n="$pairs" 'BEGIN{printf "%.6f", 1000000/n}')
    echo "[SCALE] $sample -> $scale (from $pairs spike-in pairs; $mapped_reads reads)" >&2
  else
    echo "[WARN] $sample: no spike-in pairs; using scale=1" >&2
  fi
  echo "$scale" | tee "$scale_file"
}

make_bigwigs() {
  local sample="$1"
  local bam="${sample}.filtered.bam"
  local scale; scale=$(get_spike_scale "$sample")
  local out_bw="alignment/bigwig/${sample}.spikenorm.bw"

  if [[ "$scale" != "1" && ! -f "$out_bw" ]]; then
    echo "[BIGWIG] spike-in normalized $sample"
    bamCoverage -b "$bam" -o "$out_bw" -bs "$BW_BINSIZE" \
      --scaleFactor "$scale" -p "$THREADS"
  elif [[ -f "$out_bw" ]]; then
    echo "[SKIP] spike-in bigWig exists: $out_bw"
  else
    echo "[INFO] Skipping spike-in bigWig (scale=1)."
  fi
}

make_bigwigs "$CONTROL"
for s in "${SAMPLES[@]}"; do
  make_bigwigs "$s"
done

# ==============================
# BEDGRAPH + SEACR
# ==============================
BG_PATH=""
ensure_bedgraph() {
  local sample="$1"
  local bam="${sample}.filtered.bam"
  local out="alignment/bedgraph/${sample}.spike.bedgraph"
  if [[ ! -s "$out" ]]; then
    [[ -s "$bam" ]] || { echo "[ERROR] missing $bam" >&2; return 1; }
    local scale; scale=$(get_spike_scale "$sample")
    echo "[BG] $sample: genomecov scale=$scale -> $out" >&2
    bedtools genomecov -bg -pc -ibam "$bam" \
      | awk -v s="$scale" 'BEGIN{OFS="\t"} NF==4 && $2<$3 { $4=s*$4; print }' \
      | LC_ALL=C sort -k1,1 -k2,2n > "$out"
  else
    echo "[SKIP] $out exists ($(wc -l < "$out") lines)" >&2
  fi
  BG_PATH="$out"
}

# control bedGraph
ensure_bedgraph "$CONTROL"
CTRL_BG="$BG_PATH"
echo "[INFO] Control bedGraph: $CTRL_BG" >&2

# run SEACR
for sample in "${SAMPLES[@]}"; do
  ensure_bedgraph "$sample"
  EXP_BG="$BG_PATH"

  if [[ -s "$EXP_BG" && -s "$CTRL_BG" ]]; then
    out_prefix="peakCalling/SEACR/${sample}_vsControl"
    echo "[SEACR] $sample vs $CONTROL -> $out_prefix" >&2
    bash "$SEACR_SH" "$EXP_BG" "$CTRL_BG" non stringent "$out_prefix" \
      || echo "[WARNING][$sample] SEACR vs control failed." >&2

    out_top="peakCalling/SEACR/${sample}_top0.01"
    bash "$SEACR_SH" "$EXP_BG" 0.01 non stringent "$out_top" \
      || echo "[WARNING][$sample] SEACR top0.01 failed." >&2
  else
    echo "[WARNING][$sample] Skipping SEACR (missing/empty EXP or CTRL bedGraph)." >&2
  fi
  echo "[DEBUG][$sample] Processing completed." >&2
done

echo "[DONE] IgGCAU + Irf1PCAU: Trimming + Mapping + Spike-in BigWigs + SEACR complete."
