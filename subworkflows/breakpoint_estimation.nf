#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Define genomic breakpoints for JAbBA
//

// Delly
include { BAM_DELLY } from '../subworkflows/bam_delly/main.nf'
include { BCFTOOLS_VIEW } from '../../modules/local/bcftools/view/main.nf'

// GRIDSS
include { GRIDSS_SV_CALLING } from '../subworkflows/bam_gridss/main.nf'
include { GRIDSS_SOMATIC_FILTER_STEP } from '../subworkflows/bam_gridss/main.nf'

// Genomic breakpoint estimation workflow for JAbBA
workflow BREAKPOINT_ESTIMATOR {
    take:
    ch_samples       // channel: [ val(meta), control, control_index, tumor, tumor_index ]
    fasta            
    fai
    bwa_index
    gridss_blacklist
    gridss_pon
    delly_blacklist
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
        .combine(delly_blacklist)

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

    //
    // Run GRIDSS SV calling
    //
    GRIDSS_SV_CALLING(
        ch_samples,
        bwa_index,
        fasta,
        fai,
        gridss_blacklist
    )
    versions = versions.mix(GRIDSS_SV_CALLING.out.versions)

    //
    // Run GRIDSS Somatic Filter
    //
    GRIDSS_SOMATIC_FILTER_STEP(
        GRIDSS_SV_CALLING.out.vcf,
        gridss_pon
    )
    gridss_vcf_hc   = GRIDSS_SOMATIC_FILTER_STEP.out.somatic_high_confidence
    gridss_vcf_all  = GRIDSS_SOMATIC_FILTER_STEP.out.somatic_all_vcf
    versions        = versions.mix(GRIDSS_SOMATIC_FILTER_STEP.out.versions)

    emit:
    delly_filter_bcf
    delly_vcf
    gridss_vcf_hc
    gridss_vcf_all
    versions
}