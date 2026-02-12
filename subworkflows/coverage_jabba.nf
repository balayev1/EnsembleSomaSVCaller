#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Generate GC- and mappability-corrected denoised genomic coverage files for JABBA
//

params.options = [:]

include { FRAGCOUNTER } from './modules/nf-core/fragcounter/main.nf' addParams( options: params.options )
include { REBIN_RAW_FRAGCOUNTER } from './modules/nf-core/fragcounter/main.nf' addParams( options: params.options )

// Genomic coverage acquisition workflow for JABBA
workflow COVERAGE_JABBA {
    take
        ch_samples,    // channel: [val(meta), tumor,tumor_bai, control, control_bai]
        midpoint_frag  
        windowsize_frag
        gcmapdir_frag
        minmapq_frag
        paired_frag
        exome_frag
    
    main:
    versions = Channel.empty()
    fragcounter_raw_cov_tumor   = Channel.empty()
    fragcounter_cov_tumor   = Channel.empty()
    corrected_bw_tumor      = Channel.empty()
    rebinned_raw_cov_tumor  = Channel.empty()

    // Tumor input
    ch_tumor = ch_samples.map { meta, tumor, tumor_bai, control, control_bai ->
        return [ meta, tumor, tumor_bai ]
    }

    // Control input
    ch_control = ch_samples.map { meta, tumor, tumor_bai, control, control_bai ->
        return [ meta, control, control_bai ]
    }

    // Combine inputs
    ch_fragcounter_input = ch_tumor.mix(ch_control)

    // Run fragCounter on tumor and control samples
    FRAGCOUNTER(
        ch_fragcounter_input, 
        midpoint, 
        windowsize, 
        gcmapdir, 
        minmapq, 
        fasta, 
        fasta_fai, 
        paired, 
        exome
    )

    // Initialize outputs from fragCounter
    fragcounter_raw_cov = FRAGCOUNTER.out.fragcounter_raw_cov
    fragcounter_cov = FRAGCOUNTER.out.fragcounter_cov
    corrected_bw = FRAGCOUNTER.out.corrected_bw
    versions = FRAGCOUNTER.out.versions

    // Rebin fragCounter output into 1kb windows
    REBIN_RAW_FRAGCOUNTER(fragcounter_cov, "reads", 1000)
    rebinned_raw_cov  = REBIN_RAW_FRAGCOUNTER.out.raw_fragcounter_cov_1kb

    //
    emit:
    fragcounter_raw_cov
    fragcounter_cov
    rebinned_raw_cov
    corrected_bw
    
    versions
}
