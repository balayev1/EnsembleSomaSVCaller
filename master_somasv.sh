#!/bin/bash

# master_somasv.sh
# Master script to run nf-aceseq and EnsembleSomaSVCaller pipelines on SLURM
# author: balayev1 & gemini pro & claude sonnet 4.5
# Usage: bash master_somasv.sh

# Load modules
module load nextflow
module load singularity

# TODO: SET SINGULARITY CACHE DIR
export NXF_SINGULARITY_CACHEDIR="path/to/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# TODO: DEFINE PATHS TO GITHUB REPOS
ACESEQ_REPO="path/to/nf-aceseq"
SOMASV_REPO="path/to/EnsembleSomaSVCaller"

# # ----------------------------------------------------
# # 1. Clone & Prep Repos
# # ----------------------------------------------------

# # Clone nf-aceseq
# if [ ! -d "$ACESEQ_REPO" ]; then
#     echo "Cloning nf-aceseq..."
#     git clone https://github.com/ghga-de/nf-aceseq.git
# fi
# # ALWAYS ensure bin is executable (even if repo existed previously)
# chmod +x nf-aceseq/bin/

# # Clone EnsembleSomaSVCaller
# if [ ! -d "$SOMASV_REPO" ]; then
#     echo "Cloning EnsembleSomaSVCaller..."
#     git clone https://github.com/balayev1/EnsembleSomaSVCaller.git
# fi

# ----------------------------------------------------
# 2. Set Arguments
# ----------------------------------------------------

# TODO: UPDATE & ADD ARGUMENTS AS NEEDED! CHECK REPO FOR MORE INFO!
## GENERAL ARGUMENTS
export PROFILE="singularity"

## ACESEQ-SPECIFIC ARGUMENTS 
export ACESEQ_INPUT="path/to/samplesheet.csv" 
export ACESEQ_OUTPUT="path/to/ACESEQ_out"
export BEAGLE_REF="path/to/beagle_ref"
export BEAGLE_REF_EXT="vcf"
export BEAGLE_GMAP="path/to/beagle_gmap/no_chr_in_chrom_field"
export BEAGLE_GMAP_EXT="map"
export GENOME_VER="GRCh38"
export MPILEUP_QUAL=20
export SNP_MIN_COVERAGE=15
export REFGENIE_IGNORE=true

# ----------------------------------------------------
# 3. Submit Jobs
# ----------------------------------------------------

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