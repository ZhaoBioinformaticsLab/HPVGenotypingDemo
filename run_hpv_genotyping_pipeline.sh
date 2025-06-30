#!/bin/bash
#SBATCH --job-name=hpv_genotyping
#SBATCH --output=hpv_genotyping.%j.out
#SBATCH --error=hpv_genotyping.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load required modules on Biowulf
module purge
module load python/3.10
module load bowtie2
module load samtools
module load fastp

# Optional: activate virtual environment if using one
# source /data/$USER/envs/hpv_env/bin/activate

# Run the pipeline
echo "Starting HPV genotyping pipeline..."
python ./scripts/HPV_genotyping_pipeline_demo.py
echo "Pipeline finished."
