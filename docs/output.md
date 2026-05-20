# Outputs

## Overview

Published files are organized under `--outdir`. By default, `--outdir` is derived from `--somasv_out_base`, and Nextflow execution reports are written to `pipeline_info/`.

The exact contents of each directory depend on sample success and tool behavior, but the layout below reflects the publish rules configured in `conf/modules.config`.

## Top-Level Output Structure

- `pipeline_info/`: execution report, trace, DAG, and timeline HTML or text reports
- `ascat/<sample>/`: ASCAT purity/ploidy and segmentation outputs
- `sequenza/reference/`: Sequenza reference wiggle file
- `sequenza/<sample>/`: Sequenza purity/ploidy fits, segments, and plots
- `facets/<sample>/`: FACETS pileup, VCF, and QC plots
- `purity_ploidy/<sample>/`: merged purity/ploidy candidates and consensus tables
- `hetpileups/<sample>/`: heterozygous pileup text files for JaBbA
- `delly/<sample>/call/`: raw DELLY calls
- `delly/<sample>/filter/`: filtered DELLY calls
- `delly/<sample>/view/`: DELLY VCF conversions used downstream
- `gridss/<sample>/call/`: GRIDSS calls and assembly BAM
- `gridss/<sample>/filter/`: somatic GRIDSS VCFs
- `svaba/<sample>/`: SvABA somatic and unfiltered callsets plus breakpoint text outputs
- `survivor/merge/<sample>/`: merged consensus breakpoint VCFs
- `survivor/filter/<sample>/`: filtered consensus breakpoint VCFs
- `bam_stats/normal/<sample>/`: normal-sample BAM statistics used by BRASS
- `bam_stats/tumor/<sample>/`: tumor-sample BAM statistics used by BRASS
- `brass/<sample>/`: BRASS BEDPE and VCF outputs
- `Coverages/fragCounter_tumor/<sample>/`: tumor fragCounter coverage files
- `Coverages/fragCounter_normal/<sample>/`: normal fragCounter coverage files
- `Coverages/Dryclean_tumor/<sample>/`: tumor Dryclean corrected coverage
- `Coverages/Dryclean_normal/<sample>/`: normal Dryclean corrected coverage
- `Coverages/CBS/<sample>/`: CBS segmentation inputs for JaBbA
- `jabba/<sample>/`: JaBbA graph outputs

## Key Result Files

### Purity and ploidy

Important files published during purity/ploidy estimation include:

- ASCAT: `*purityploidy.txt`, `*segments.txt`, `*metrics.txt`, `*png`
- Sequenza: `*_alternative_solutions.txt`, `*_segments.txt`, `*_genome_view.pdf`, `*_model_fit.pdf`
- FACETS: `*.vcf.gz`, `*.csv.gz`, `*.cnv.png`, `*.cov.pdf`, `*.spider.pdf`
- Consensus merge: `*.purity_ploidy_candidates.tsv`, `*.purity_ploidy_consensus.tsv`

The consensus table is the file consumed downstream to select JaBbA purity and ploidy unless explicit overrides are supplied.

Purity/ploidy consensus is estimated per sample from candidate `(purity, ploidy)` points. The candidate sources are ASCAT, FACETS, Sequenza, and ACEseq by default; when `--skip_aceseq true` is set, ACEseq is omitted and the consensus is based on ASCAT, FACETS, and Sequenza only.

The merge step clusters candidate points in two-dimensional purity/ploidy space with DBSCAN. By default, points are grouped when they are within `eps = 0.1`, and a cluster must contain at least `min_samples = 2` points. A cluster is only eligible for consensus if it is supported by at least two distinct base callers. Multiple ACEseq rows count as multiple candidate points, but they still represent one base caller, so ACEseq alone cannot form a valid consensus cluster.

If multiple eligible clusters are found, the selected cluster is ranked by more supporting callers, then more supporting points, then tighter clustering around its centroid. The reported consensus purity and ploidy are the medians of the purity and ploidy values in that selected cluster.

At least two expected callers must provide usable purity/ploidy values. With ACEseq enabled, this means at least two of ASCAT, FACETS, Sequenza, and ACEseq; with `--skip_aceseq true`, at least two of ASCAT, FACETS, and Sequenza. If too many callers are missing or no eligible cluster is found, `*.purity_ploidy_consensus.tsv` reports `status = no_consensus` with `purity = NA` and `ploidy = NA`.

### Breakpoint calling and consensus SVs

Breakpoint evidence is generated from multiple callers and harmonized before JaBbA:

- DELLY: raw and filtered VCF or BCF outputs
- GRIDSS: raw call VCFs, indexes, and assembly BAM
- SvABA: somatic SV VCFs and breakpoint support text files
- BRASS: BEDPE and VCF outputs
- SURVIVOR: merged and filtered consensus breakpoint VCFs

The filtered SURVIVOR VCF is the principal breakpoint input passed to JaBbA.

Before JaBbA, consensus breakpoint records are annotated with a support-derived `TIER` value. Structural variants supported by two or more breakpoint callers are assigned `TIER=1`; variants supported by only one caller are assigned `TIER=2`. This tiering is based on SURVIVOR caller support after merging and filtering.

### Coverage products

Coverage processing produces:

- fragCounter raw and corrected coverage RDS files
- optional bigWig tracks from fragCounter
- Dryclean corrected coverage RDS files
- CBS `*.cov.rds`, `*.seg.rds`, and `*.nseg.rds`

The CBS segmentation outputs are passed directly into JaBbA together with the filtered consensus junctions.

### JaBbA outputs

The JaBbA publish directory contains the main graph products:

- `*.jabba.simple.rds`
- `*.jabba.simple.gg.rds`
- `*.jabba.simple.vcf`
- `*.jabba.raw.rds`
- `*opt.report.rds`
- `*.jabba.seg`
- `*karyograph.rds`

These are the final graph reconstruction outputs for each sample.

## Versions And Provenance

Most published tool directories also include a `versions.yml` file capturing software versions for that process. Combined with the `pipeline_info/` reports, these files provide a lightweight provenance trail for each run.
