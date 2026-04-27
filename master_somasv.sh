#!/bin/bash

set -euo pipefail

module load singularity

usage() {
    cat <<'USAGE'
Usage:
  bash master_somasv.sh \
    --somasv-repo /path/to/EnsembleSomaSVCaller \
    --aceseq-repo /path/to/nf-aceseq \
    --somasv-out-base /path/to/SomaticSV_outs \
    --samplesheet /path/to/samplesheet.csv \
    [--run-id RUN_ID] \
    [--gurobi-path /path/to/gurobi.lic]

Required:
  --somasv-repo
  --aceseq-repo
  --somasv-out-base
  --samplesheet

Optional:
  --run-id
  --gurobi-path

Legacy positional form is still accepted:
  bash master_somasv.sh /path/to/samplesheet.csv [run_id]
USAGE
}

SOMASV_REPO="${SOMASV_REPO:-}"
ACESEQ_REPO="${ACESEQ_REPO:-}"
SOMASV_OUT_BASE="${SOMASV_OUT_BASE:-}"
SAMPLESHEET="${SAMPLESHEET:-}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
GUROBI_PATH="${GUROBI_PATH:-}"

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --somasv-repo)
            SOMASV_REPO="$2"
            shift 2
            ;;
        --aceseq-repo)
            ACESEQ_REPO="$2"
            shift 2
            ;;
        --somasv-out-base)
            SOMASV_OUT_BASE="$2"
            shift 2
            ;;
        --samplesheet)
            SAMPLESHEET="$2"
            shift 2
            ;;
        --run-id)
            RUN_ID="$2"
            shift 2
            ;;
        --gurobi-path)
            GUROBI_PATH="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        --*)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

if [[ -z "${SAMPLESHEET}" && ${#POSITIONAL[@]} -ge 1 ]]; then
    SAMPLESHEET="${POSITIONAL[0]}"
fi
if [[ ${#POSITIONAL[@]} -ge 2 ]]; then
    RUN_ID="${POSITIONAL[1]}"
fi

if [[ -z "${SOMASV_REPO}" || -z "${ACESEQ_REPO}" || -z "${SOMASV_OUT_BASE}" || -z "${SAMPLESHEET}" ]]; then
    usage >&2
    exit 1
fi

if [[ ! -d "${SOMASV_REPO}" ]]; then
    echo "SOMASV_REPO not found: ${SOMASV_REPO}" >&2
    exit 1
fi
if [[ ! -f "${SOMASV_REPO}/main.nf" ]]; then
    echo "SOMASV_REPO does not contain main.nf: ${SOMASV_REPO}" >&2
    exit 1
fi
if [[ ! -d "${ACESEQ_REPO}" ]]; then
    echo "ACESEQ_REPO not found: ${ACESEQ_REPO}" >&2
    exit 1
fi
if [[ ! -f "${ACESEQ_REPO}/main.nf" ]]; then
    echo "ACESEQ_REPO does not contain main.nf: ${ACESEQ_REPO}" >&2
    exit 1
fi
if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "Samplesheet not found: ${SAMPLESHEET}" >&2
    exit 1
fi

if [[ -z "${GUROBI_PATH}" ]]; then
    GUROBI_PATH="${SOMASV_REPO%/}/assets/gurobi.lic"
fi
if [[ ! -f "${GUROBI_PATH}" ]]; then
    echo "Gurobi license not found: ${GUROBI_PATH}" >&2
    exit 1
fi

export SOMASV_REPO
export ACESEQ_REPO
export SOMASV_OUT_BASE
export GUROBI_PATH

ACESEQ_OUT_BASE="${ACESEQ_OUT_BASE:-${SOMASV_OUT_BASE%/}/ACESEQ_out}"
ACESEQ_RUN_OUTDIR="${ACESEQ_OUT_BASE%/}/${RUN_ID}"
SOMASV_RUN_OUTDIR="${SOMASV_OUT_BASE%/}/SOMASV_out/${RUN_ID}"
HANDOFF_DIR="${ACESEQ_RUN_OUTDIR}/handoff"
LAUNCH_ROOT="${SOMASV_OUT_BASE%/}/launches/${RUN_ID}"

mkdir -p \
    "${SOMASV_REPO}/logs" \
    "${ACESEQ_RUN_OUTDIR}" \
    "${SOMASV_RUN_OUTDIR}" \
    "${HANDOFF_DIR}" \
    "${LAUNCH_ROOT}/samplesheets"



cd "${SOMASV_REPO}"

echo "Run ID: ${RUN_ID}"
echo "Samplesheet: ${SAMPLESHEET}"
echo "ACEseq outdir: ${ACESEQ_RUN_OUTDIR}"
echo "SomaSV outdir: ${SOMASV_RUN_OUTDIR}"
echo "Gurobi license: ${GUROBI_PATH}"
echo "Launch root: ${LAUNCH_ROOT}"
echo

{
    IFS= read -r header || {
        echo "Failed to read header from samplesheet: ${SAMPLESHEET}" >&2
        exit 1
    }

    while IFS= read -r line; do
        [[ -z "${line}" ]] && continue

        sample="$(printf '%s\n' "${line}" | cut -d',' -f1)"
        sample_sheet="${LAUNCH_ROOT}/samplesheets/${sample}.csv"
        aceseq_manifest="${HANDOFF_DIR}/${sample}/aceseq_manifest.tsv"
        aceseq_launch_dir="${LAUNCH_ROOT}/aceseq/${sample}"
        somasv_launch_dir="${LAUNCH_ROOT}/somasv/${sample}"
        aceseq_work_dir="${LAUNCH_ROOT}/work/aceseq/${sample}"
        somasv_work_dir="${LAUNCH_ROOT}/work/somasv/${sample}"

        mkdir -p \
            "$(dirname "${aceseq_manifest}")" \
            "${aceseq_work_dir}" \
            "${somasv_work_dir}"

        {
            printf '%s\n' "${header}"
            printf '%s\n' "${line}"
        } > "${sample_sheet}"

        ACESEQ_JOB=$(
            sbatch --parsable \
                --export=ALL,SOMASV_REPO="${SOMASV_REPO}",ACESEQ_REPO="${ACESEQ_REPO}",SOMASV_OUT_BASE="${SOMASV_OUT_BASE}",PIPELINE_INPUT="${sample_sheet}",PIPELINE_SAMPLE_ID="${sample}",PIPELINE_OUTDIR="${ACESEQ_RUN_OUTDIR}",PIPELINE_WORKDIR="${aceseq_work_dir}",PIPELINE_LAUNCH_DIR="${aceseq_launch_dir}",ACESEQ_MANIFEST="${aceseq_manifest}" \
                "${SOMASV_REPO}/slurm/nf_aceseq.sbatch"
        )

        SOMASV_JOB=$(
            sbatch --parsable \
                --dependency=afterok:${ACESEQ_JOB} \
                --export=ALL,SOMASV_REPO="${SOMASV_REPO}",SOMASV_OUT_BASE="${SOMASV_OUT_BASE}",PIPELINE_INPUT="${sample_sheet}",PIPELINE_SAMPLE_ID="${sample}",PIPELINE_OUTDIR="${SOMASV_RUN_OUTDIR}",PIPELINE_WORKDIR="${somasv_work_dir}",PIPELINE_LAUNCH_DIR="${somasv_launch_dir}",ACESEQ_MANIFEST="${aceseq_manifest}",PIPELINE_GUROBI_PATH="${GUROBI_PATH}" \
                "${SOMASV_REPO}/slurm/nf_somasv.sbatch"
        )

        printf 'Sample %-20s ACEseq job %-12s SomaSV job %-12s Manifest %s\n' \
            "${sample}" "${ACESEQ_JOB}" "${SOMASV_JOB}" "${aceseq_manifest}"
    done
} < "${SAMPLESHEET}"
