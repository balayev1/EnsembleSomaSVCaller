#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Generate GC- and mappability-corrected denoised genomic coverage files for JAbBA
//

params.options = [:]

// fragCounter
include { BAM_FRAGCOUNTER as TUMOR_FRAGCOUNTER         } from './subworkflows/bam_fragCounter/main'
include { BAM_FRAGCOUNTER as NORMAL_FRAGCOUNTER        } from './subworkflows/bam_fragCounter/main'

// dryclean
include { COV_DRYCLEAN as TUMOR_DRYCLEAN               } from './subworkflows/cov_dryclean/main'
include { COV_DRYCLEAN as NORMAL_DRYCLEAN              } from './subworkflows/cov_dryclean/main'

// Genomic coverage acquisition workflow for JAbBA
workflow COVERAGE_JABBA {
    take
        ch_samples,    // channel: [ val(meta), [control],[control_index],[ tumor], [tumor_index]]
        midpoint_frag  
        windowsize_frag
        gcmapdir_frag
        minmapq_frag
        paired_frag
        exome_frag
        pon_dryclean
        center_dryclean
        cbs_dryclean
        cnsignif_dryclean
        wholeGenome_dryclean
        blacklist_dryclean
        blacklist_path_dryclean
        germline_filter_dryclean
        germline_file_dryclean
        field_dryclean
        build_dryclean
    
    main:
    versions = Channel.empty()

    // Control input
    ch_control = ch_samples.map { meta, control, control_index, tumor, tumor_index ->
        return [ meta, control, control_index ]
    }

    // Tumor input
    ch_tumor = ch_samples.map { meta, control, control_index, tumor, tumor_index ->
        return [ meta, tumor, tumor_index ]
    }

    // Run fragCounter for tumor and normal samples to GC- and mappability-correct coverage for JAbBA
    NORMAL_FRAGCOUNTER(ch_control, midpoint_frag, windowsize_frag, gcmapdir_frag, minmapq_frag, paired_frag, exome_frag)
    //normal_frag_cov = Channel.empty().mix(NORMAL_FRAGCOUNTER.out.fragcounter_cov)
    normal_frag_cov = Channel.empty().mix(NORMAL_FRAGCOUNTER.out.rebinned_raw_cov)

    TUMOR_FRAGCOUNTER(ch_tumor, midpoint_frag, windowsize_frag, gcmapdir_frag, minmapq_frag, paired_frag, exome_frag)
    //tumor_frag_cov = Channel.empty().mix(TUMOR_FRAGCOUNTER.out.fragcounter_cov)
    tumor_frag_cov = Channel.empty().mix(TUMOR_FRAGCOUNTER.out.rebinned_raw_cov)

    // Acquire fragCounter version
    versions = versions.mix(NORMAL_FRAGCOUNTER.out.versions)

    // Run dryclean for tumor and normal samples to denoise genomic coverage for JAbBA
    TUMOR_DRYCLEAN(tumor_frag_cov, pon_dryclean, center_dryclean,
                    cbs_dryclean, cnsignif_dryclean, wholeGenome_dryclean,
                    blacklist_dryclean, blacklist_path_dryclean,
                    germline_filter_dryclean, germline_file_dryclean, 
                    field_dryclean, build_dryclean)

    tumor_dryclean_cov = Channel.empty().mix(TUMOR_DRYCLEAN.out.dryclean_cov)

    NORMAL_DRYCLEAN(normal_frag_cov, pon_dryclean, center_dryclean,
            cbs_dryclean, cnsignif_dryclean, wholeGenome_dryclean,
            blacklist_dryclean, blacklist_path_dryclean,
            germline_filter_dryclean, germline_file_dryclean,
            field_dryclean, build_dryclean)

    normal_dryclean_cov = Channel.empty().mix(NORMAL_DRYCLEAN.out.dryclean_cov)

    // Acquire dryclean version
    versions = versions.mix(TUMOR_DRYCLEAN.out.versions) 

    //
    emit:
    tumor_dryclean_cov  = tumor_dryclean_cov
    normal_dryclean_cov = normal_dryclean_cov
    versions            = versions
}
