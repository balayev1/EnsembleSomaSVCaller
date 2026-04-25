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

set -euo pipefail

unset JAVA_HOME JAVA_CMD
module purge
module load java/openjdk-21.0.2
module load nextflow/25.10.4
module load singularity

PIPELINE_INPUT="${PIPELINE_INPUT:?PIPELINE_INPUT is required}"
PIPELINE_SAMPLE_ID="${PIPELINE_SAMPLE_ID:?PIPELINE_SAMPLE_ID is required}"
PIPELINE_OUTDIR="${PIPELINE_OUTDIR:?PIPELINE_OUTDIR is required}"
PIPELINE_WORKDIR="${PIPELINE_WORKDIR:?PIPELINE_WORKDIR is required}"
PIPELINE_LAUNCH_DIR="${PIPELINE_LAUNCH_DIR:?PIPELINE_LAUNCH_DIR is required}"
ACESEQ_MANIFEST="${ACESEQ_MANIFEST:?ACESEQ_MANIFEST is required}"
SOMASV_REPO="${SOMASV_REPO:?SOMASV_REPO is required}"
PIPELINE_GUROBI_PATH="${PIPELINE_GUROBI_PATH:-}"

mkdir -p "${PIPELINE_LAUNCH_DIR}" "${PIPELINE_WORKDIR}"
cd "${PIPELINE_LAUNCH_DIR}"

echo "SomaSV Pipeline Started at $(date)"

python3 "${SOMASV_REPO}/bin/validate_aceseq_manifest.py" \
    --samplesheet "${PIPELINE_INPUT}" \
    --manifest "${ACESEQ_MANIFEST}"

nextflow_cmd=(
    nextflow run "${SOMASV_REPO}/main.nf"
    -profile singularity
    --input "${PIPELINE_INPUT}"
    --outdir "${PIPELINE_OUTDIR}"
    --aceseq_manifest "${ACESEQ_MANIFEST}"
    -work-dir "${PIPELINE_WORKDIR}"
    -resume
)

if [[ -n "${PIPELINE_GUROBI_PATH}" ]]; then
    nextflow_cmd+=( --gurobi_license_path "${PIPELINE_GUROBI_PATH}" )
fi

"${nextflow_cmd[@]}"

echo "SomaSV Pipeline Finished at $(date)"
