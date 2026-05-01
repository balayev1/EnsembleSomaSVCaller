# Documentation

This directory contains the user-facing documentation for EnsembleSomaSVCaller.

## Contents

- [usage.md](usage.md): how to prepare inputs, launch the pipeline, and override references or runtime settings
- [output.md](output.md): published output layout and key result files produced for each sample

## Pipeline Summary

EnsembleSomaSVCaller processes matched tumor-normal BAM files and ACEseq-derived auxiliary inputs to produce somatic structural variant calls and JaBbA genome graph reconstructions. At a high level, the pipeline performs the following stages:

1. Validate the main samplesheet and ACEseq manifest.
2. Normalize dbSNP contig naming for downstream compatibility.
3. Estimate purity and ploidy with ASCAT, Sequenza, FACETS, and ACEseq-derived inputs.
4. Build GC- and mappability-corrected coverage tracks for JaBbA.
5. Generate heterozygous pileups for JaBbA.
6. Call breakpoints with DELLY, GRIDSS, SvABA, and BRASS.
7. Merge and filter breakpoint evidence with SURVIVOR.
8. Segment corrected coverage with CBS.
9. Run JaBbA using consensus breakpoints, segmented coverage, het pileups, and purity/ploidy inputs.

## Notes

- Parameter validation and samplesheet validation are enabled through `nf-schema`.
- Most reference resources are resolved from `params.genomes` using the configured `--genome` key.
- This repository is a custom research pipeline, not a full nf-core template pipeline, so the included lint configuration intentionally disables a small set of template-specific checks.
