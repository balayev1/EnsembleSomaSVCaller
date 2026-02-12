#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Generate GC- and mappability-corrected denoised genomic coverage files for JABBA
//

params.options = [:]

// fragCounter
include { BAM_FRAGCOUNTER as TUMOR_FRAGCOUNTER         } from '../subworkflows/bam_fragCounter/main'
include { BAM_FRAGCOUNTER as NORMAL_FRAGCOUNTER        } from '../subworkflows/bam_fragCounter/main'

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
    fragcounter_raw_cov   = Channel.empty()
    fragcounter_cov_tumor   = Channel.empty()
    corrected_bw_tumor      = Channel.empty()
    rebinned_raw_cov_tumor  = Channel.empty()

    // Control input
    ch_control = ch_samples.map { meta, tumor, tumor_bai, control, control_bai ->
        return [ meta, control, control_bai ]
    }

    // Tumor input
    ch_tumor = ch_samples.map { meta, tumor, tumor_bai, control, control_bai ->
        return [ meta, tumor, tumor_bai ]
    }

    // Run fragCounter for tumor and normal samples
    NORMAL_FRAGCOUNTER(ch_control, midpoint_frag, windowsize_frag, gcmapdir_frag, minmapq_frag, paired_frag, exome_frag)
    //normal_frag_cov = Channel.empty().mix(NORMAL_FRAGCOUNTER.out.fragcounter_cov)
    normal_frag_cov = Channel.empty().mix(NORMAL_FRAGCOUNTER.out.rebinned_raw_cov)

    TUMOR_FRAGCOUNTER(ch_tumor, midpoint_frag, windowsize_frag, gcmapdir_frag, minmapq_frag, paired_frag, exome_frag)
    //tumor_frag_cov = Channel.empty().mix(TUMOR_FRAGCOUNTER.out.fragcounter_cov)
    tumor_frag_cov = Channel.empty().mix(TUMOR_FRAGCOUNTER.out.rebinned_raw_cov)

    // Acquire fragCounter version
    versions = versions.mix(NORMAL_FRAGCOUNTER.out.versions)


}
