#!/bin/bash

# master_somasv.sh
# Master script to run nf-aceseq and EnsembleSomaSVCaller pipelines on SLURM
# author: balayev1 & gemini pro & claude sonnet 4.5
# Usage: bash master_somasv.sh

# Load modules
module load nextflow
module load singularity

# TODO: DEFINE PATHS TO GITHUB REPOS
ACESEQ_REPO="path/to/nf-aceseq"
SOMASV_REPO="path/to/EnsembleSomaSVCaller"

# Create log directories for BOTH pipelines
mkdir -p "${ACESEQ_REPO}/logs"
mkdir -p "${SOMASV_REPO}/logs"

# --- JOB 1: nf-aceseq ---
cd "$ACESEQ_REPO" || exit
ACESEQ_JOB=$(sbatch --parsable "${SOMASV_REPO}/slurm/nf_aceseq.sbatch")
echo "Submitted ACEseq job: $ACESEQ_JOB"

# --- JOB 2: SomaSV (Dependent) ---
cd "$SOMASV_REPO" || exit
SOMASV_JOB=$(sbatch --parsable --dependency=afterok:${ACESEQ_JOB} "${SOMASV_REPO}/slurm/nf_somasv.sbatch")
echo "Submitted SomaSV job: $SOMASV_JOB (Waiting for $ACESEQ_JOB)"