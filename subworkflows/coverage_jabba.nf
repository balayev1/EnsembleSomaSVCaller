#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Generate GC- and mappability-corrected denoised genomic coverage files for JABBA
//

include { BAM_FRAGCOUNTER as TUMOR_FRAGCOUNTER }  from './bam_fragCounter/main.nf'
include { BAM_FRAGCOUNTER as NORMAL_FRAGCOUNTER } from './bam_fragCounter/main.nf'
include { COV_DRYCLEAN as TUMOR_DRYCLEAN }        from './cov_dryclean/main.nf'
include { COV_DRYCLEAN as NORMAL_DRYCLEAN }       from './cov_dryclean/main.nf'

workflow COVERAGE_JABBA {
    take:
        ch_samples
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

        // INPUT_PREP emits [meta, control, control_index, tumor, tumor_index]
        ch_control = ch_samples.map { meta, control, control_bai, tumor, tumor_bai ->
            [meta, control, control_bai]
        }

        ch_tumor = ch_samples.map { meta, control, control_bai, tumor, tumor_bai ->
            [meta, tumor, tumor_bai]
        }

        NORMAL_FRAGCOUNTER(
            ch_control,
            midpoint_frag,
            windowsize_frag,
            gcmapdir_frag,
            minmapq_frag,
            paired_frag,
            exome_frag
        )
        normal_frag_cov = NORMAL_FRAGCOUNTER.out.rebinned_raw_cov

        TUMOR_FRAGCOUNTER(
            ch_tumor,
            midpoint_frag,
            windowsize_frag,
            gcmapdir_frag,
            minmapq_frag,
            paired_frag,
            exome_frag
        )
        tumor_frag_cov = TUMOR_FRAGCOUNTER.out.rebinned_raw_cov

        versions = versions
            .mix(NORMAL_FRAGCOUNTER.out.versions)
            .mix(TUMOR_FRAGCOUNTER.out.versions)

        TUMOR_DRYCLEAN(
            tumor_frag_cov,
            pon_dryclean,
            center_dryclean,
            cbs_dryclean,
            cnsignif_dryclean,
            wholeGenome_dryclean,
            blacklist_dryclean,
            blacklist_path_dryclean,
            germline_filter_dryclean,
            germline_file_dryclean,
            field_dryclean,
            build_dryclean
        )
        tumor_dryclean_cov = TUMOR_DRYCLEAN.out.dryclean_cov

        NORMAL_DRYCLEAN(
            normal_frag_cov,
            pon_dryclean,
            center_dryclean,
            cbs_dryclean,
            cnsignif_dryclean,
            wholeGenome_dryclean,
            blacklist_dryclean,
            blacklist_path_dryclean,
            germline_filter_dryclean,
            germline_file_dryclean,
            field_dryclean,
            build_dryclean
        )
        normal_dryclean_cov = NORMAL_DRYCLEAN.out.dryclean_cov

        versions = versions
            .mix(TUMOR_DRYCLEAN.out.versions)
            .mix(NORMAL_DRYCLEAN.out.versions)

    emit:
        tumor_dryclean_cov
        normal_dryclean_cov
        versions
}
