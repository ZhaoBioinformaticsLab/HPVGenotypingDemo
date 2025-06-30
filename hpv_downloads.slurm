#!/bin/bash
#SBATCH --job-name=hpv_sra_download
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --gres=lscratch:20
#SBATCH --time=04:00:00
#SBATCH --output=hpv_download.%j.out
#SBATCH --error=hpv_download.%j.err

# Load SRA Toolkit
module load sratoolkit

# Print node and environment info
echo "[START] Running on $(hostname)"
echo "[START] Working directory: $(pwd)"
echo "[START] Job ID: $SLURM_JOB_ID"

# Use current directory
cd . || exit 1

# Check script existence
if [[ ! -x download_hpv_amplicons.sh ]]; then
  echo "[ERROR] Script download_hpv_amplicons.sh not found or not executable"
  exit 1
fi

# Run the download script
bash download_hpv_amplicons.sh

echo "[DONE] SLURM job completed"
