#!/bin/bash
#SBATCH --job-name=SomaSV_Master
#SBATCH --partition=asvnode1,msibigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nf_somasv_%j.out
#SBATCH --error=logs/nf_somasv_%j.err
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

echo "SomaSV Pipeline Started at $(date)"

# Run the pipeline
nextflow run main.nf \
    -profile singularity \
    -resume

echo "SomaSV Pipeline Finished at $(date)"