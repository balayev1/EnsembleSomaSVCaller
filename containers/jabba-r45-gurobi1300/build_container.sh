#!/bin/bash
#SBATCH --job-name=genomic_analysis
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=128gb
#SBATCH --partition=msibigmem,msismall,msilarge,agsmall
#SBATCH --mail-type=ALL
#SBATCH --mail-user=balay011@umn.edu
#SBATCH -e /scratch.global/venteicher_30050/balay011/EnsembleSomaSVCaller/containers/jabba-r45-gurobi1300/%x_%j.err
#SBATCH -o /scratch.global/venteicher_30050/balay011/EnsembleSomaSVCaller/containers/jabba-r45-gurobi1300/%x_%j.out

module load apptainer
export APPTAINER_TMPDIR=/scratch.global/venteicher_30050/balay011/EnsembleSomaSVCaller/containers/jabba-r45-gurobi1300/apptainer_tmp
export APPTAINER_CACHEDIR=/scratch.global/venteicher_30050/balay011/singularity_cache

mkdir -p "${APPTAINER_TMPDIR}"
cd /scratch.global/venteicher_30050/balay011/EnsembleSomaSVCaller/containers/jabba-r45-gurobi1300

apptainer build --fakeroot \
	jabba-r45-gurobi1300_0.2.0.sif \
	apptainer.def
