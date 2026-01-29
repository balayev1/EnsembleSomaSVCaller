#!/bin/bash
#SBATCH --job-name=ASCAT_Master
#SBATCH --partition=asvnode1,msibigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=ascat_master_%j.out
#SBATCH --error=ascat_master_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

# Load modules
module load nextflow
module load singularity
module load java/openjdk-17.0.2

# Run Nextflow
nextflow run main.nf \
    -profile singularity \
    -resume \
    -c nextflow.config \
    -with-report report.html \
    -with-timeline timeline.html

echo "ASCAT Pipeline Finished at $(date)"