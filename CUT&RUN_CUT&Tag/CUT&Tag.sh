#!/bin/bash
#SBATCH --job-name CUTTag_pipeline
#SBATCH --error CUTTag_pipeline-error.e%j
#SBATCH --output CUTTag_pipeline-out.o%j
#SBATCH --cpus-per-task 64
#SBATCH --mem=128000
#SBATCH --partition shared-cpu
#SBATCH --time 11:00:00

###############################################################################
#                             PIPELINE SETTINGS
###############################################################################
set -uo pipefail

# Species: "mouse" or "human"
SPECIES="human"

# Prefix of control IgG sample (must match FASTQ basename before first underscore)
CONTROL_IGG="RabbitIgGCAU"

# Toggle duplicate removal: "yes" or "no"
USE_DEDUP="yes"

# Toggle E. coli spike-in normalization: "yes" or "no"
USE_SPIKEIN="yes"

# Path to E. coli Bowtie2 index (prefix only, no .bt2)
SPIKEIN_INDEX="/home/users/j/javed/scratch/Fastqs_050525/tools/Ecoli_genome/genome"

# FASTQ directory and OUTPUT directory
FASTQ_DIR="/srv/beegfs/scratch/users/j/javed/Fastqs_050525/allfastqs/CT_Org"
OUTPUT_DIR="/home/users/j/javed/scratch/Fastqs_050525/allfastqs/CT_Org/outs_dupstandalone"

# SEACR path
SEACR="/home/users/j/javed/scratch/Fastqs_050525/tools/SEACR_1.3.sh"

# Threads
THREADS=${SLURM_CPUS_PER_TASK:-8}

###############################################################################
#             REFERENCES, EFFECTIVE GENOME SIZE & BLACKLIST FILES
###############################################################################
if [[ "$SPECIES" == "mouse" ]]; then
  genomeIndex="/home/users/j/javed/scratch/Fastqs_050525/tools/mm10arc"
  chromSize="/home/users/j/javed/scratch/Fastqs_050525/tools/mm10arcchrom.txt"
  effectiveGenomeSize=2652783500
  blacklist="/home/users/j/javed/scratch/Fastqs_050525/tools/mm10-blacklist.v2.bed"
else
  genomeIndex="/home/users/j/javed/scratch/Fastqs_050525/tools/GRCh38_noalt_as/GRCh38_noalt_as"
  chromSize="/home/users/j/javed/scratch/Fastqs_050525/tools/hg38.chrom.sizes"
  effectiveGenomeSize=2913022398
  blacklist="/home/users/j/javed/scratch/Fastqs_050525/tools/hg38-blacklist.v2.bed"
fi

###############################################################################
#                               LOAD MODULES
###############################################################################
module load picard/2.21.1-Java-11

# Configure Picard invocation via EBROOTPICARD
: "${EBROOTPICARD:?Picard module didnt set EBROOTPICARD}"
[[ -f "$EBROOTPICARD/picard.jar" ]] || { echo "[ERROR] $EBROOTPICARD/picard.jar not found"; exit 1; }

MEM_MB=${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_TASK:-128000}}
JAVA_XMX_GB=$(( MEM_MB / 1024 - 8 ))
if (( JAVA_XMX_GB < 4 )); then JAVA_XMX_GB=4; fi
PICARD_CMD=(java -Xmx${JAVA_XMX_GB}g -jar "$EBROOTPICARD/picard.jar")

###############################################################################
#                       CREATE OUTPUT DIRECTORIES
###############################################################################
mkdir -p \
  "$OUTPUT_DIR/alignment/sam" \
  "$OUTPUT_DIR/alignment/bowtie2_summary" \
  "$OUTPUT_DIR/alignment/removeDuplicate" \
  "$OUTPUT_DIR/alignment/bam" \
  "$OUTPUT_DIR/alignment/bed" \
  "$OUTPUT_DIR/alignment/bedgraph" \
  "$OUTPUT_DIR/alignment/bigwig" \
  "$OUTPUT_DIR/alignment/filtered_bam" \
  "$OUTPUT_DIR/alignment/spikein" \
  "$OUTPUT_DIR/tmp" \
  "$OUTPUT_DIR/peakCalling/SEACR"

echo "========== PIPELINE STARTED at $(date) =========="

###############################################################################
# 1) DISCOVER SAMPLE PREFIXES
###############################################################################
echo "[DEBUG] Discovering sample names..."
mapfile -t allSamples < <(
  ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null \
    | xargs -n1 basename \
    | cut -d'_' -f1 \
    | sort -u
)
if (( ${#allSamples[@]} == 0 )); then
  echo "[ERROR] No FASTQs found in $FASTQ_DIR"
  exit 1
fi
echo "[DEBUG] Found ${#allSamples[@]} samples: ${allSamples[*]}"

###############################################################################
# 2) ALIGNMENT FUNCTION (genome + optional spike-in)
###############################################################################
align_reads() {
  local sample=$1
  echo "[DEBUG][$sample] Starting Bowtie2 alignment..."

  local R1s=( "$FASTQ_DIR"/${sample}_*_R1_001.fastq.gz )
  local R2s=( "$FASTQ_DIR"/${sample}_*_R2_001.fastq.gz )

  if (( ${#R1s[@]} == 0 || ${#R2s[@]} == 0 )); then
    echo "[ERROR][$sample] No FASTQs found. Skipping..."
    return 1
  fi

  local R1list; R1list=$(printf ",%s" "${R1s[@]}"); R1list=${R1list:1}
  local R2list; R2list=$(printf ",%s" "${R2s[@]}"); R2list=${R2list:1}

  local outSam="$OUTPUT_DIR/alignment/sam/${sample}.sam"
  local log="$OUTPUT_DIR/alignment/bowtie2_summary/${sample}.log"

  if ! bowtie2 -p "$THREADS" -x "$genomeIndex" -1 "$R1list" -2 "$R2list" -S "$outSam" 2> "$log"; then
    echo "[ERROR][$sample] Bowtie2 failed. Check $log"
    return 1
  fi
  echo "[DEBUG][$sample] Bowtie2 alignment completed."

  # SPIKE-IN ALIGNMENT
  if [[ "$USE_SPIKEIN" == "yes" ]]; then
    echo "[DEBUG][$sample] Starting spike-in alignment..."
    local spikeSam="$OUTPUT_DIR/alignment/spikein/${sample}_spikein.sam"
    local spikeLog="$OUTPUT_DIR/alignment/spikein/${sample}_spikein.log"

    if bowtie2 --end-to-end --very-sensitive \
               --no-overlap --no-dovetail --no-mixed --no-discordant \
               --no-unal -I 10 -X 700 \
               -x "$SPIKEIN_INDEX" \
               -1 "$R1list" -2 "$R2list" \
               -S "$spikeSam" \
               -p 8 2> "$spikeLog"; then

      local spike_aligns
      spike_aligns=$(samtools view -F 0x04 "$spikeSam" | wc -l || echo 0)
      local spike_pairs=$(( spike_aligns / 2 ))
      echo "$spike_pairs" > "$OUTPUT_DIR/alignment/spikein/${sample}_spikein.seqDepth"

      if (( spike_pairs == 0 )); then
        echo "[WARNING][$sample] No spike-in reads. Will fall back to RPGC scaling."
      else
        echo "[DEBUG][$sample] Spike-in mapped pairs: $spike_pairs"
      fi
    else
      echo "[WARNING][$sample] Spike-in alignment failed. Falling back to RPGC."
      echo 0 > "$OUTPUT_DIR/alignment/spikein/${sample}_spikein.seqDepth"
    fi
  fi

  return 0
}

###############################################################################
# 3) SPIKE-IN SCALE FACTOR
###############################################################################
compute_spikein_scale() {
  local sample=$1
  local sfile="$OUTPUT_DIR/alignment/spikein/${sample}_spikein.seqDepth"
  local spikeDepth=0
  [[ -f "$sfile" ]] && spikeDepth=$(<"$sfile")

  local minfile="$OUTPUT_DIR/alignment/spikein/_GLOBAL_MIN.seqDepth"
  local min_spike=0
  if [[ -f "$minfile" ]]; then
    min_spike=$(<"$minfile")
  else
    local tmpmin
    tmpmin=$(ls "$OUTPUT_DIR/alignment/spikein/"*.seqDepth 2>/dev/null | xargs -n1 cat | sort -n | head -1)
    min_spike=${tmpmin:-0}
  fi

  if [[ -z "${spikeDepth:-}" || -z "${min_spike:-}" || "$spikeDepth" -eq 0 || "$min_spike" -eq 0 ]]; then
    echo 1
  else
    echo "$(echo "$min_spike / $spikeDepth" | bc -l)"
  fi
}

###############################################################################
# 4) PROCESS ONE SAMPLE (Picard MarkDuplicates → filter → coverage → SEACR)
###############################################################################
process_sample() {
  local sample=$1
  local sam="$OUTPUT_DIR/alignment/sam/${sample}.sam"
  [[ -f "$sam" ]] || { echo "[ERROR][$sample] SAM missing. Skipping..."; return 1; }

  echo "[DEBUG][$sample] Processing..."

  # Convert SAM -> coordinate-sorted BAM
  local possorted="$OUTPUT_DIR/alignment/bam/${sample}.positionsort.bam"
  if ! samtools view -@ "$THREADS" -bS "$sam" | samtools sort -@ "$THREADS" -o "$possorted" -; then
    echo "[ERROR][$sample] Convert + coord-sort failed."
    return 1
  fi

  # Picard MarkDuplicates via java -jar $EBROOTPICARD/picard.jar
  local dedup="$OUTPUT_DIR/alignment/removeDuplicate/${sample}.rmdup.bam"
  local metrics="$OUTPUT_DIR/alignment/removeDuplicate/${sample}.markdup.metrics.txt"
  local finalbam

  if [[ "$USE_DEDUP" == "yes" ]]; then
    if ! "${PICARD_CMD[@]}" MarkDuplicates \
        I="$possorted" \
        O="$dedup" \
        M="$metrics" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        ASSUME_SORTED=true \
        TMP_DIR="$OUTPUT_DIR/tmp"; then
      echo "[ERROR][$sample] Picard MarkDuplicates failed."
      return 1
    fi
    finalbam="$dedup"
  else
    finalbam="$possorted"
  fi

  # Check mapped reads pre-blacklist
  local mapped_pre
  mapped_pre=$(samtools view -c -F 0x4 "$finalbam" || echo 0)
  if (( mapped_pre == 0 )); then
    echo "[ERROR][$sample] No mapped reads after (de)dup. Skipping downstream."
    return 1
  fi

  # Blacklist filtering (re-sort to ensure coordinate order)
  local filtered_bam="$OUTPUT_DIR/alignment/filtered_bam/${sample}.filtered.bam"
  if ! bedtools intersect -v -abam "$finalbam" -b "$blacklist" \
      | samtools sort -@ "$THREADS" -o "$filtered_bam" -; then
    echo "[ERROR][$sample] Blacklist filtering failed."
    return 1
  fi
  samtools index -@ "$THREADS" "$filtered_bam" || true

  local mapped_filt
  mapped_filt=$(samtools view -c -F 0x4 "$filtered_bam" || echo 0)
  if (( mapped_filt == 0 )); then
    echo "[ERROR][$sample] No mapped reads after blacklist filter. Skipping downstream."
    return 1
  fi

  # BED + bedGraph (sorted) for SEACR
  local filtbed="$OUTPUT_DIR/alignment/bed/${sample}_filtered.bed"
  if ! bedtools bamtobed -i "$filtered_bam" > "$filtbed"; then
    echo "[ERROR][$sample] bamtobed failed."
    return 1
  fi
  if [[ ! -s "$filtbed" ]]; then
    echo "[ERROR][$sample] Empty BED. Skipping coverage/SEACR."
    return 1
  fi

  local bg scale seqDepth
  if [[ "$USE_SPIKEIN" == "yes" ]]; then
    scale=$(compute_spikein_scale "$sample")
    bg="$OUTPUT_DIR/alignment/bedgraph/${sample}_scaled.bedgraph"
    if ! bedtools genomecov -bg -scale "$scale" -i "$filtbed" -g "$chromSize" > "$bg"; then
      echo "[ERROR][$sample] genomecov (spike-in scaled) failed."
      return 1
    fi
  else
    local aligns
    aligns=$(samtools view -c -F 0x4 "$filtered_bam" || echo 0)
    seqDepth=$(( aligns / 2 ))
    if (( seqDepth == 0 )); then
      echo "[ERROR][$sample] seqDepth=0 (RPGC). Skipping."
      return 1
    fi
    bg="$OUTPUT_DIR/alignment/bedgraph/${sample}_rpgc.bedgraph"
    local factor
    factor=$(echo "1 / ($seqDepth / 10000000)" | bc -l)
    if ! bedtools genomecov -bg -scale "$factor" -i "$filtbed" -g "$chromSize" > "$bg"; then
      echo "[ERROR][$sample] genomecov (RPGC) failed."
      return 1
    fi
  fi

  # Ensure bedGraph coordinate order for SEACR
  LC_ALL=C sort -k1,1 -k2,2n -o "$bg" "$bg"

  # BigWigs
  if [[ "$USE_SPIKEIN" == "yes" ]]; then
    bamCoverage -b "$filtered_bam" -o "$OUTPUT_DIR/alignment/bigwig/${sample}.bw" --binSize 1 -p "$THREADS" --scaleFactor "$scale" |
| {
      echo "[WARNING][$sample] bamCoverage (scaled) failed."; }
    bamCoverage -b "$filtered_bam" -o "$OUTPUT_DIR/alignment/bigwig/${sample}_raw.bw" --binSize 1 -p "$THREADS" || {
      echo "[WARNING][$sample] bamCoverage (raw) failed."; }
  else
    bamCoverage -b "$filtered_bam" -o "$OUTPUT_DIR/alignment/bigwig/${sample}.bw" --binSize 1 -p "$THREADS" --normalizeUsing RPGC --
effectiveGenomeSize "$effectiveGenomeSize" || {
      echo "[WARNING][$sample] bamCoverage (RPGC) failed."; }
    bamCoverage -b "$filtered_bam" -o "$OUTPUT_DIR/alignment/bigwig/${sample}_raw.bw" --binSize 1 -p "$THREADS" || {
      echo "[WARNING][$sample] bamCoverage (raw) failed."; }
  fi
  echo "[DEBUG][$sample] BigWigs generated."

  # SEACR (skip if inputs missing/empty)
  local ctrl_bg
  if [[ "$USE_SPIKEIN" == "yes" ]]; then
    ctrl_bg="$OUTPUT_DIR/alignment/bedgraph/${CONTROL_IGG}_scaled.bedgraph"
  else
    ctrl_bg="$OUTPUT_DIR/alignment/bedgraph/${CONTROL_IGG}_rpgc.bedgraph"
  fi

  if [[ -s "$bg" && -s "$ctrl_bg" ]]; then
    bash "$SEACR" "$bg" "$ctrl_bg" non stringent "$OUTPUT_DIR/peakCalling/SEACR/${sample}_vsControl.peaks" || \
      echo "[WARNING][$sample] SEACR vs control failed."
    bash "$SEACR" "$bg" 0.01 non stringent "$OUTPUT_DIR/peakCalling/SEACR/${sample}_top0.01.peaks" || \
      echo "[WARNING][$sample] SEACR top0.01 failed."
  else
    echo "[WARNING][$sample] Skipping SEACR (experiment or control bedGraph empty/missing)."
  fi

  echo "[DEBUG][$sample] Processing completed."
  return 0
}

###############################################################################
# 5) RUN: ALIGN ALL SAMPLES FIRST → THEN PROCESS ALL
###############################################################################
# Align all (so spike-in depths exist for all before scaling)
for smp in "${allSamples[@]}"; do
  if ! align_reads "$smp"; then
    echo "[WARNING][$smp] Alignment failed. Will skip processing."
  fi
done

# Global min spike-in depth (for stable scaling)
if [[ "$USE_SPIKEIN" == "yes" ]]; then
  GLOBAL_MIN_SPIKE=$(ls "$OUTPUT_DIR/alignment/spikein/"*.seqDepth 2>/dev/null | xargs -n1 cat | sort -n | head -1)
  : > "$OUTPUT_DIR/alignment/spikein/_GLOBAL_MIN.seqDepth"
  echo "${GLOBAL_MIN_SPIKE:-0}" > "$OUTPUT_DIR/alignment/spikein/_GLOBAL_MIN.seqDepth"
  echo "[DEBUG] Global min spike-in pairs: ${GLOBAL_MIN_SPIKE:-0}"
fi

# Process control first (so ctrl bedGraph exists for SEACR)
if [[ " ${allSamples[*]} " == *" $CONTROL_IGG "* ]]; then
  process_sample "$CONTROL_IGG" || echo "[WARNING][$CONTROL_IGG] Processing failed."
fi

# Process all samples
for smp in "${allSamples[@]}"; do
  [[ "$smp" == "$CONTROL_IGG" ]] && continue
  process_sample "$smp" || echo "[WARNING][$smp] Processing failed."
done

echo "========== PIPELINE FINISHED at $(date) =========="
