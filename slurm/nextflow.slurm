#!/bin/bash
#SBATCH --job-name=SomaSV_Master
#SBATCH --partition=asvnode1,msibigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nf_somasv_%j.out
#SBATCH --error=logs/nf_somasv_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

# Make logs directory
mkdir -p logs

# Load modules
module load nextflow
module load singularity
module load java/openjdk-17.0.2

# Singularity cachedir
export NXF_SINGULARITY_CACHEDIR="/scratch.global/balay011/singularity_cache"
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Run the pipeline
nextflow run main.nf \
    -profile singularity \
    -resume

echo "SomaSV Pipeline Finished at $(date)"