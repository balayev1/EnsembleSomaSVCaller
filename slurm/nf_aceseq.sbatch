#!/bin/bash
#SBATCH --job-name=nf_aceseq
#SBATCH --partition=asvnode1,msibigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nf_aceseq_%j.out
#SBATCH --error=logs/nf_aceseq_%j.err
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

echo "nf-aceseq started at $(date)"

nextflow run main.nf \
    --input "$ACESEQ_INPUT" \
    --outdir "$ACESEQ_OUTPUT" \
    -profile "$PROFILE" \
    --mpileup_qual "$MPILEUP_QUAL" \
    --snp_min_coverage "$SNP_MIN_COVERAGE" \
    --genome "$GENOME_VER" \
    -resume

echo "nf-aceseq finished at $(date)"