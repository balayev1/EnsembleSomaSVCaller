#!/bin/bash
#SBATCH --job-name=nf_aceseq
#SBATCH --partition=asvnode1,msibigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --output=logs/nf_aceseq_%j.out
#SBATCH --error=logs/nf_aceseq_%j.err
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
ACESEQ_REPO="${ACESEQ_REPO:?ACESEQ_REPO is required}"
SOMASV_REPO="${SOMASV_REPO:?SOMASV_REPO is required}"
ACESEQ_EXTRA_CONFIG="${ACESEQ_EXTRA_CONFIG:-}"

extra_config_args=()
if [[ -n "${ACESEQ_EXTRA_CONFIG}" ]]; then
    extra_config_args=( -c "${ACESEQ_EXTRA_CONFIG}" )
fi

mkdir -p "${PIPELINE_LAUNCH_DIR}" "${PIPELINE_WORKDIR}" "$(dirname "${ACESEQ_MANIFEST}")"
cd "${PIPELINE_LAUNCH_DIR}"

echo "nf-aceseq started at $(date)"

nextflow run "${ACESEQ_REPO}/main.nf" \
    -profile singularity \
    -c "${SOMASV_REPO}/conf/aceseq_msi.config" \
    "${extra_config_args[@]}" \
    --input "${PIPELINE_INPUT}" \
    --outdir "${PIPELINE_OUTDIR}" \
    --skip_multiqc true \
    -work-dir "${PIPELINE_WORKDIR}" \
    -resume

python3 "${SOMASV_REPO}/bin/build_aceseq_manifest.py" \
    --sample-id "${PIPELINE_SAMPLE_ID}" \
    --aceseq-root "${PIPELINE_OUTDIR}" \
    --output "${ACESEQ_MANIFEST}"

python3 "${SOMASV_REPO}/bin/validate_aceseq_manifest.py" \
    --samplesheet "${PIPELINE_INPUT}" \
    --manifest "${ACESEQ_MANIFEST}"

echo "nf-aceseq finished at $(date)"
