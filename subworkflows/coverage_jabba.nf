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


    //
    // Run fragCounter on tumor and control samples
    //

}
