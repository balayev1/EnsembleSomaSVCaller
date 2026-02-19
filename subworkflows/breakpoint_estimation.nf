#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Define genomic breakpoints for JAbBA
//

// Delly
include { BAM_DELLY } from '../subworkflows/bam_delly/main.nf'
include { BCFTOOLS_VIEW } from '../../modules/local/bcftools/view/main.nf'

// Genomic breakpoint estimation workflow for JAbBA
workflow BREAKPOINT_ESTIMATOR {
    take:
    ch_samples       // channel: [ val(meta), control, control_index, tumor, tumor_index ]
    fasta            // channel: [ val(meta2), path(fasta) ]
    fai              // channel: [ val(meta3), path(fai) ]
    exclude_bed
    delly_mode
    delly_min_svqual
    delly_altaf
    delly_min_svsize
    delly_max_svsize
    delly_ratio_geno
    delly_pass
    delly_tags
    delly_cov
    delly_ctrl_contamination
    delly_geno_qual
    delly_rddel
    delly_rddup

    main:
    versions = Channel.empty()

    // Prepare reference inputs for DELLY
    ch_fasta = fasta.map { it -> [ [id:'fasta'], it ] }
    ch_fai   = fai.map   { it -> [ [id:'fasta_fai'], it ] }

    // Transform ch_samples into DELLY input format
    delly_input = ch_samples
        .map { meta, control, control_index, tumor, tumor_index ->
            [ meta, [control, tumor], [control_index, tumor_index] ]
        }
        .combine(exclude_bed)

    //
    // Run DELLY
    //
    BAM_DELLY (
        delly_input,
        ch_fasta,
        ch_fai,
        delly_mode,
        delly_min_svqual,
        delly_altaf,
        delly_min_svsize,
        delly_max_svsize,
        delly_ratio_geno,
        delly_pass,
        delly_tags,
        delly_cov,
        delly_ctrl_contamination,
        delly_geno_qual,
        delly_rddel,
        delly_rddup
    )
    delly_filter_bcf        = BAM_DELLY.out.delly_filter_bcf
    versions                = versions.mix(BAM_DELLY.out.versions)

    //
    // Run bcftools view to convert DELLY BCF to VCF for merging
    //
    bcftools_input = delly_filter_bcf.map { meta, bcf -> [ meta, bcf, [] ] }
    BCFTOOLS_VIEW (
        bcftools_input,
        [],
        [],
        []
    )
    
    delly_vcf = BCFTOOLS_VIEW.out.vcf
    versions  = versions.mix(BCFTOOLS_VIEW.out.versions_bcftools)

    emit:
    delly_filter_bcf
    delly_vcf
    versions = versions
}