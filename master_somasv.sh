#!/bin/bash

# Load modules
module load nextflow
module load singularity

# Singularity cachedir
export NXF_SINGULARITY_CACHEDIR="/scratch.global/balay011/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# Define github repo paths (Define these early so we can use them)
BASE_DIR="$PWD"
ACESEQ_REPO="${BASE_DIR}/nf-aceseq"
SOMASV_REPO="${BASE_DIR}/EnsembleSomaSVCaller"

# ----------------------------------------------------
# 1. Clone & Prep Repos
# ----------------------------------------------------

# Clone nf-aceseq
if [ ! -d "$ACESEQ_REPO" ]; then
    echo "Cloning nf-aceseq..."
    git clone https://github.com/ghga-de/nf-aceseq.git
fi
# ALWAYS ensure bin is executable (even if repo existed previously)
chmod +x "${ACESEQ_REPO}/bin/"*

# Clone EnsembleSomaSVCaller
if [ ! -d "$SOMASV_REPO" ]; then
    echo "Cloning EnsembleSomaSVCaller..."
    git clone https://github.com/balayev1/EnsembleSomaSVCaller.git
fi

# ----------------------------------------------------
# 2. Set Arguments
# ----------------------------------------------------

# TODO: UPDATE & ADD AS NEEDED!
export ACESEQ_INPUT="${BASE_DIR}/samplesheet.csv" 
export ACESEQ_OUTPUT="${BASE_DIR}/ACESEQ_out"

export PROFILE="singularity"
export MPILEUP_QUAL=20
export SNP_MIN_COVERAGE=15
export GENOME_VER="GRCh38"

# ----------------------------------------------------
# 3. Submit Jobs
# ----------------------------------------------------

# Create log directories for BOTH pipelines
mkdir -p "${ACESEQ_REPO}/logs"
mkdir -p "${SOMASV_REPO}/logs"

# --- JOB 1: nf-aceseq ---
cd "$ACESEQ_REPO" || exit
ACESEQ_JOB=$(sbatch --parsable "${SOMASV_REPO}/nf_aceseq.sbatch")
echo "Submitted ACEseq job: $ACESEQ_JOB"

# --- JOB 2: SomaSV (Dependent) ---
cd "$SOMASV_REPO" || exit
SOMASV_JOB=$(sbatch --parsable --dependency=afterok:${ACESEQ_JOB} "${SOMASV_REPO}/nf_somasv.sbatch")
echo "Submitted SomaSV job: $SOMASV_JOB (Waiting for $ACESEQ_JOB)"