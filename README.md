# EnsembleSomaSVCaller

EnsembleSomaSVCaller is a Nextflow DSL2 pipeline for somatic structural variant calling with JaBbA-based genome graph reconstruction from matched tumor-normal whole-genome sequencing data.

This repository is designed to work downstream of [`nf-aceseq`](https://github.com/ghga-de/nf-aceseq), which produces the ACEseq purity/ploidy and segmentation outputs that EnsembleSomaSVCaller consumes during purity/ploidy reconciliation and JaBbA preparation.

If needed, the pipeline can also run in a reduced `--skip_aceseq true` mode. In that mode, ACEseq is omitted entirely and purity/ploidy consensus is computed from ASCAT, FACETS, and Sequenza only.

## End-To-End Workflow

The default intended flow is:

1. Prepare one shared tumor-normal samplesheet.
2. Run `nf-aceseq` on the same samples.
3. Build an `aceseq_manifest.tsv` that points to the required ACEseq outputs.
4. Validate the manifest against the samplesheet.
5. Run EnsembleSomaSVCaller with both the samplesheet and the ACEseq manifest.

Alternative flow:

1. Prepare one shared tumor-normal samplesheet.
2. Launch EnsembleSomaSVCaller with `--skip_aceseq true`.
3. The pipeline skips ACEseq manifest validation and merges purity/ploidy from the remaining callers only.

If you are on MSI/SLURM and want that chaining done for you, this repo also includes [`master_somasv.sh`](master_somasv.sh), [`slurm/nf_aceseq.sbatch`](slurm/nf_aceseq.sbatch), and [`slurm/nf_somasv.sbatch`](slurm/nf_somasv.sbatch).

## Repositories

Clone both pipelines:

```bash
git clone https://github.com/balayev1/EnsembleSomaSVCaller.git
git clone https://github.com/ghga-de/nf-aceseq.git
```

## Requirements

- Nextflow `25.10.2` or newer
- Java 17 or newer
- Singularity or Apptainer
- Matched tumor-normal BAMs and BAM indexes
- Access to all reference resources required by both pipelines
- Recommended for production JaBbA runs: a valid Gurobi license plus a reachable Gurobi installation that can be exposed as `GUROBI_HOME`
- A CPLEX installation exposed as `CPLEXDIR` only if you intentionally want JaBbA's non-Gurobi `Rcplex` fallback

Important local note:

- [`conf/aceseq_msi.config`](conf/aceseq_msi.config) contains MSI-specific reference and output defaults. Review and update those paths before using the bundled ACEseq launch helper on another system.
- [`nextflow.config`](nextflow.config) contains the default EnsembleSomaSVCaller reference bundle and runtime defaults. For portable usage, prefer explicit CLI parameters or a custom `-c` config instead of editing project defaults repeatedly.
- The current project default in [`nextflow.config`](nextflow.config) is `gurobi_jabba = false`. If you want the recommended Gurobi path, pass `--gurobi_jabba TRUE` explicitly or override it in your own config.

## Shared Samplesheet

Both pipelines use the same CSV samplesheet structure:

- `sample`: unique sample or case identifier used consistently across the samplesheet, ACEseq manifest, and published output directories
- `sex`: biological sex label used by the copy-number workflows, typically `male` or `female`
- `tumor`: path to the tumor BAM
- `tumor_index`: path to the tumor BAM index, usually `*.bai`
- `control`: path to the matched normal or control BAM
- `control_index`: path to the matched normal or control BAM index, usually `*.bai`

Optional:

- `sv`: optional legacy metadata field stored as `meta.sv_type`; it is not required by the downstream EnsembleSomaSVCaller steps currently implemented in this repository

Example:

```csv
sample,sex,sv,tumor,tumor_index,control,control_index
CB822,male,,/data/CB822_tumor.bam,/data/CB822_tumor.bam.bai,/data/CB822_normal.bam,/data/CB822_normal.bam.bai
```

## Step 1: Run nf-aceseq

Upstream `nf-aceseq` documentation describes a standard Nextflow launch using `--input`, `--outdir`, and a container profile. In this repo, the bundled MSI helper runs `nf-aceseq` with the local ACEseq config in [`conf/aceseq_msi.config`](conf/aceseq_msi.config), which is a good starting point if your environment is similar.

Optional skip-ACEseq route: if ACEseq outputs are unavailable or you want to run the reduced caller ensemble, skip Step 1 through Step 3 and launch EnsembleSomaSVCaller with `--skip_aceseq true` in Step 4. In this mode the pipeline does not require `--aceseq_manifest`; purity/ploidy consensus is computed from ASCAT, FACETS, and Sequenza only.

Example direct run:

```bash
nextflow run /path/to/nf-aceseq/main.nf \
  -profile singularity \
  -c /path/to/EnsembleSomaSVCaller/conf/aceseq_msi.config \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/ACESEQ_out
```

If your `nf-aceseq` checkout does not preserve executable bits in `bin/`, upstream also notes that you may need:

```bash
chmod +x /path/to/nf-aceseq/bin/*
```

## Step 2: Build The ACEseq Manifest

EnsembleSomaSVCaller does not ingest the full `nf-aceseq` directory blindly. Instead, it expects an `aceseq_manifest.tsv` with one row per sample and the following columns:

- `sample_id`
- `aceseq_dir`
- `cnv_tab`
- `snp_tab`
- `pscbs_data`
- `sex_file`
- `clustered_and_pruned_and_normal`
- `ploidy_purity_2d`
- `breakpoints2`
- `known_segments`

This repo includes [`bin/build_aceseq_manifest.py`](bin/build_aceseq_manifest.py), which builds a one-sample manifest row from an ACEseq output directory.

Example for one sample:

```bash
python3 bin/build_aceseq_manifest.py \
  --sample-id CB822 \
  --aceseq-root /path/to/ACESEQ_out \
  --output /path/to/CB822.aceseq_manifest.tsv
```

The script expects the following `nf-aceseq` sample outputs to exist:

- `{sample}.cnv.tab.gz`
- `{sample}.snp.tab.gz`
- `{sample}_pscbs_data.txt.gz`
- `{sample}_sex.txt`
- `{sample}_clustered_and_pruned_and_normal.txt`
- `{sample}_ploidy_purity_2D.txt`
- `{sample}_breakpoints2.txt`
- `{sample}.knownSegments.txt`

For a cohort-level EnsembleSomaSVCaller run, combine the per-sample rows into one manifest with a single header. For example:

```bash
{
  head -n 1 /path/to/CB822.aceseq_manifest.tsv
  tail -n +2 /path/to/*.aceseq_manifest.tsv
} > /path/to/aceseq_manifest.tsv
```

## Step 3: Validate The ACEseq Manifest

Before launching EnsembleSomaSVCaller, validate that the manifest covers every sample in the main samplesheet and that all referenced ACEseq artifacts exist:

```bash
python3 bin/validate_aceseq_manifest.py \
  --samplesheet /path/to/samplesheet.csv \
  --manifest /path/to/aceseq_manifest.tsv
```

Skip this step if you plan to run with `--skip_aceseq true`.

## Step 4: Run EnsembleSomaSVCaller

Launch EnsembleSomaSVCaller with the same samplesheet plus the validated ACEseq manifest:

```bash
nextflow run main.nf \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results
```

To run without ACEseq:

```bash
nextflow run main.nf \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --skip_aceseq true \
  --outdir /path/to/results
```

Recommended JaBbA solver path: Gurobi

For most real samples, run JaBbA with Gurobi enabled. The JaBbA module expects:

- `GUROBI_HOME` to point to the Gurobi installation root, either from your shell environment or from a custom Nextflow config that sets `params.gurobi_home`
- `--gurobi_license_path` to point to a valid `gurobi.lic`

Example direct Nextflow launch:

```bash
export GUROBI_HOME=/opt/gurobi1300/linux64

nextflow run main.nf \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results \
  --gurobi_jabba TRUE \
  --gurobi_license_path /path/to/gurobi.lic
```

Fallback path: JaBbA without Gurobi (`Rcplex`/CPLEX)

If you set `--gurobi_jabba FALSE`, the pipeline passes `--gurobi FALSE` to JaBbA and uses the `Rcplex`/CPLEX code path instead. That path requires `CPLEXDIR` to be defined, either in your shell environment or through a custom Nextflow config that sets `params.cplex_dir`.

Example direct Nextflow launch:

```bash
export CPLEXDIR=/opt/cplex2211

nextflow run main.nf \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results \
  --gurobi_jabba FALSE
```

Important note:

- the bundled JaBbA image is `oras://quay.io/balay011/jabba-r45-gurobi1300:0.2.0`
- that same image may be used for both JaBbA with Gurobi and JaBbA with the `Rcplex`/CPLEX fallback path
- in practice, the fallback path may still fail on realistic samples if CPLEX is running in Community Edition or other size-restricted mode
- if you see `CPLEX Error 1016: Community Edition. Problem size limits exceeded`, switch to the Gurobi path for production runs or provide a full non-size-restricted CPLEX setup


## One-Command MSI/SLURM Launch

If you want this repository to submit `nf-aceseq` and EnsembleSomaSVCaller in sequence, use [`master_somasv.sh`](master_somasv.sh):

```bash
bash master_somasv.sh \
  --somasv-repo /path/to/EnsembleSomaSVCaller \
  --aceseq-repo /path/to/nf-aceseq \
  --somasv-out-base /path/to/SomaticSV_outs \
  --samplesheet /path/to/samplesheet.csv \
  --gurobi-jabba TRUE \
  --gurobi-home /opt/gurobi1300/linux64 \
  --gurobi-path /path/to/gurobi.lic
```

To submit SomaSV directly without ACEseq through the wrapper:

```bash
bash master_somasv.sh \
  --somasv-repo /path/to/EnsembleSomaSVCaller \
  --somasv-out-base /path/to/SomaticSV_outs \
  --samplesheet /path/to/samplesheet.csv \
  --skip-aceseq TRUE \
  --gurobi-jabba TRUE \
  --gurobi-home /opt/gurobi1300/linux64 \
  --gurobi-path /path/to/gurobi.lic
```

If you intentionally want the `Rcplex`/CPLEX fallback through the wrapper, disable Gurobi and provide `CPLEXDIR`:

```bash
bash master_somasv.sh \
  --somasv-repo /path/to/EnsembleSomaSVCaller \
  --aceseq-repo /path/to/nf-aceseq \
  --somasv-out-base /path/to/SomaticSV_outs \
  --samplesheet /path/to/samplesheet.csv \
  --gurobi-jabba FALSE \
  --cplex-dir /opt/cplex2211
```

What each wrapper option means:

- `--somasv-repo`: path to this `EnsembleSomaSVCaller` repository checkout. The wrapper uses it to find `main.nf`, helper scripts, and the SLURM submission scripts.
- `--aceseq-repo`: path to the `nf-aceseq` repository checkout. Required for the default chained mode; not needed with `--skip-aceseq TRUE`.
- `--skip-aceseq`: when `TRUE`, skip the upstream ACEseq submission and launch EnsembleSomaSVCaller directly without an ACEseq manifest.
- `--somasv-out-base`: base output directory for the whole chained run. The wrapper creates run-specific `ACESEQ_out/`, `SOMASV_out/`, `launches/`, and handoff directories underneath it.
- `--samplesheet`: CSV file listing the tumor-normal samples to process. The wrapper splits this into one-sample CSVs before submitting jobs.
- `--gurobi-home`: path to the Gurobi installation root. The wrapper exports it as `GUROBI_HOME` so the JaBbA task can find the solver libraries.
- `--cplex-dir`: path to the CPLEX installation root. The wrapper exports it as `CPLEXDIR` for fallback `Rcplex` runs.
- `--gurobi-jabba`: whether JaBbA should be launched with Gurobi enabled in EnsembleSomaSVCaller. The wrapper default is `FALSE`, so pass `TRUE` explicitly for the recommended production path.
- `--gurobi-path`: path to the Gurobi license file. This is only needed when `--gurobi-jabba TRUE`.

What this wrapper does:

- splits the input CSV into one sample per launch
- submits `nf-aceseq` first with [`slurm/nf_aceseq.sbatch`](slurm/nf_aceseq.sbatch) in the default mode
- builds and validates `handoff/<sample>/aceseq_manifest.tsv` in the default mode
- submits EnsembleSomaSVCaller with an `afterok` dependency using [`slurm/nf_somasv.sbatch`](slurm/nf_somasv.sbatch), or directly when `--skip-aceseq TRUE`
- forwards `--gurobi-home`, `--cplex-dir`, `--gurobi-jabba`, and `--gurobi-path` into the SLURM environment used by the downstream Nextflow JaBbA task

Typical output layout under `--somasv-out-base`:

- `ACESEQ_out/<run_id>/`
- `ACESEQ_out/<run_id>/handoff/<sample>/aceseq_manifest.tsv`
- `SOMASV_out/<run_id>/`
- `launches/<run_id>/samplesheets/`
- `launches/<run_id>/work/aceseq/<sample>/`
- `launches/<run_id>/work/somasv/<sample>/`

## Input Validation

This pipeline validates:

- top-level parameters with [`nextflow_schema.json`](nextflow_schema.json)
- the main samplesheet with [`assets/schema_input.json`](assets/schema_input.json)
- the ACEseq manifest with [`assets/schema_aceseq_manifest.json`](assets/schema_aceseq_manifest.json)

Print the generated parameter help with:

```bash
nextflow run main.nf --help
```

## Testing

The repository includes a smoke-test profile with bundled dummy inputs:

```bash
nextflow run main.nf -profile test -stub-run
```

This validates workflow wiring and schema integration, not real biological execution.

## Additional Documentation

- Pipeline overview: [docs/README.md](docs/README.md)
- Usage guide: [docs/usage.md](docs/usage.md)
- Output guide: [docs/output.md](docs/output.md)
- Citation guidance: [CITATIONS.md](CITATIONS.md)
