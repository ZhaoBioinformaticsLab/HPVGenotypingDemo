#!/bin/bash

module load sratoolkit
THREADS=8
MAIN_OUTDIR="."

# List of SRR accessions to download
SRR_LIST=(
  SRR1186008
  SRR1186009
  SRR1186010
  SRR1186011
  SRR1186012
  SRR1186013
  SRR1186014
  SRR1186015
  SRR1186016
)

for SRR in "${SRR_LIST[@]}"; do
  echo "[INFO] Processing $SRR..."

  OUTDIR="${MAIN_OUTDIR}/${SRR}"
  mkdir -p "$OUTDIR"

  # Force download directly into OUTDIR with explicit filename
  SRA_FILE="${OUTDIR}/${SRR}.sra"
  prefetch "$SRR" --output-file "$SRA_FILE"

  # Check if SRA file was downloaded successfully
  if [[ ! -s "$SRA_FILE" ]]; then
    echo "[ERROR] Failed to download $SRR"
    continue
  fi

  echo "[INFO] Converting $SRR.sra to FASTQ..."

  # Convert .sra to FASTQ
  fasterq-dump "$SRA_FILE" \
    --outdir "$OUTDIR" \
    --threads "$THREADS" \
    --split-files \
    --skip-technical \
    --include-technical NO \
    --temp /lscratch/$SLURM_JOB_ID

  echo "✔ Completed $SRR → FASTQ in $OUTDIR"
done

echo "[DONE] All downloads and conversions finished successfully.]"
