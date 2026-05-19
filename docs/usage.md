# Usage

## Requirements

- Nextflow `25.10.2` or newer
- Java 17 or newer
- A container runtime. The pipeline currently ships a `singularity` profile and is typically run with Singularity or Apptainer on HPC.
- Access to the reference bundle pointed to by `--igenomes_base`, or a custom config that overrides `params.genomes`

## Required Inputs

The pipeline always requires the main samplesheet. By default it also requires an ACEseq manifest, unless you launch with `--skip_aceseq true`.

### Main samplesheet

Pass with `--input`. This must be a CSV with the following required columns:

- `sample`
- `sex`
- `control`
- `control_index`
- `tumor`
- `tumor_index`

Optional column:

- `sv`

Example:

```csv
sample,sex,sv,tumor,tumor_index,control,control_index
CB822,male,,/data/CB822_tumor.bam,/data/CB822_tumor.bam.bai,/data/CB822_normal.bam,/data/CB822_normal.bam.bai
```

### ACEseq manifest

Pass with `--aceseq_manifest`. This must be a TSV with the following required columns. It is required unless `--skip_aceseq true` is set:

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

Each `sample_id` must match a `sample` value in the main samplesheet.

Example:

```tsv
sample_id	aceseq_dir	cnv_tab	snp_tab	pscbs_data	sex_file	clustered_and_pruned_and_normal	ploidy_purity_2d	breakpoints2	known_segments
CB822	/path/to/CB822	/path/to/CB822/cnv.tab	/path/to/CB822/snp.tab	/path/to/CB822/pscbs.RData	/path/to/CB822/sex.txt	/path/to/CB822/clustered_and_pruned_and_normal.tsv	/path/to/CB822/ploidy_purity_2d.tsv	/path/to/CB822/breakpoints2.tsv	/path/to/CB822/known_segments.tsv
```

## Running The Pipeline

### From a local clone

```bash
nextflow run main.nf \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results
```

Skip-ACEseq mode:

```bash
nextflow run main.nf \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --skip_aceseq true \
  --outdir /path/to/results
```

### From GitHub

```bash
nextflow run balayev1/EnsembleSomaSVCaller \
  -r main \
  -profile singularity \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results
```

## Parameter Validation

The pipeline validates user parameters at startup with `nf-schema`. This includes:

- top-level parameter validation from `nextflow_schema.json`
- main samplesheet validation from `assets/schema_input.json`
- ACEseq manifest validation from `assets/schema_aceseq_manifest.json` when `skip_aceseq` is not enabled

To print the generated help text:

```bash
nextflow run main.nf --help
```

## Test Profile

The pipeline includes a built-in smoke-test profile with bundled dummy manifests and placeholder reference resources:

```bash
nextflow run main.nf -profile test -stub-run
```

Notes:

- `-stub-run` is still required. The `test` profile supplies test inputs and local execution settings, but stub mode is what prevents the workflow from running the real bioinformatics tools.
- The bundled `test` profile enables Singularity so process stubs can still report tool versions from their containers.

## Common Runtime Options

- `--genome`: genome key used to resolve references from `params.genomes`
- `--igenomes_base`: base directory containing the reference bundle
- `--skip_aceseq`: skip ACEseq manifest usage and merge purity/ploidy from ASCAT, FACETS, and Sequenza only
- `--somasv_out_base`: default base directory used to derive `--outdir`
- `--publish_dir_mode`: publish mode for final outputs, default `copy`
- `--purity_jabba`: optional explicit JaBbA purity override
- `--ploidy_jabba`: optional explicit JaBbA ploidy override
- `--pp_method_jabba`: JaBbA purity/ploidy method, currently `sequenza` or `ppgrid`

## Reference Overrides

The default references are loaded from `conf/igenomes.config` through `params.genomes`. For portable runs, prefer one of these approaches instead of editing the bundled config:

### Override only the reference base

If your bundle follows the same directory structure:

```bash
nextflow run main.nf \
  -profile singularity \
  --igenomes_base /path/to/references \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results
```

### Supply a custom config

If your paths differ from the bundled layout, define a custom `params.genomes` map in a separate config and launch with `-c`:

```bash
nextflow run main.nf \
  -profile singularity \
  -c custom_genomes.config \
  --input /path/to/samplesheet.csv \
  --aceseq_manifest /path/to/aceseq_manifest.tsv \
  --outdir /path/to/results
```

## Operational Notes

- The repository currently contains site-specific defaults in `nextflow.config`. For production runs, explicitly pass `--input`, `--outdir`, and either `--aceseq_manifest` or `--skip_aceseq true`.
- The pipeline expects matched tumor-normal BAMs with valid index files.
- JaBbA purity and ploidy overrides are expected to be single values or `NA`, not candidate ranges.
