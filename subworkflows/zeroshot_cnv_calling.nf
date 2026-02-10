#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Calling CNVs using a zero-shot approach with FOUR tools ASCAT, SEQUENZA, FACETS, ACESEQ
//

params.options = [:]

include { ASCAT }                 from './modules/local/ascat' addParams( options: params.options )
include { FACETS }                from './modules/local/facets' addParams( options: params.options )
include { SEQUENZAUTILS_GCWIGGLE } from './modules/local/sequenza_utils' addParams( options: params.options )
include { SEQUENZA_RUN }          from './modules/local/sequenza_run' addParams( options: params.options )
include { CHECK_ACESEQ_DIR }      from '../modules/local/check_aceseq'

// Zero-shot CNV Calling Workflow
workflow ZERO_SHOT_CNV_CALL {
    take
        ch_samples,    // channel: [val(meta), tumor,tumor_bai, control, control_bai]
        refgen,         // channel: [path(fasta), path(fasta_fai)]
        allele_resource,    // channel: [path(allele_res)]
        loci_resource,      // channel: [path(loci_res)]
        ascat_gc_file,       // channel: [path(gc_file)]
        ascat_rt_file,       // channel: [path(rt_file)]
        ascat_bed_file,       // channel: [path(bed)]
        facets_snp, // channel: [path(facets_snp_vcf), path(facets_snp_vcf_index)]
        facets_targets, // channel: [path(facets_targets_bed)]
        facets_annotation // channel: [path(facets_annotation_bed)]

    main:
    versions = Channel.empty()

    //
    // Run ASCAT
    //
    ASCAT (
        ch_samples,
        allele_resource,
        loci_resource,
        ascat_bed_file,
        refgen[0],
        ascat_gc_file,
        ascat_rt_file
    )
    versions = versions.mix(ASCAT.out.versions)

    //
    // Run Sequenzautils to generate GC Wiggle Reference
    //
    ch_gc_wiggle_input = Channel.of( [ [id:'reference_genome'], refgen[0] ] )
    SEQUENZAUTILS_GCWIGGLE(ch_gc_wiggle_input)
    ch_wiggle_file = SEQUENZAUTILS_GCWIGGLE.out.wig.map { meta, wig -> wig }
    versions = versions.mix(SEQUENZAUTILS_GCWIGGLE.out.versions)

    // 
    // Run SEQUENZA
    //
    SEQUENZA_RUN ( 
        ch_samples, 
        refgen[0], 
        ch_wiggle_file 
    )
    versions = versions.mix(SEQUENZA_RUN.out.versions)

    // 
    // Run FACETS
    //
    FACETS (
        ch_samples,
        facets_snp[0],
        facets_snp[1],
        facets_targets,
        facets_annotation
    )
    versions = versions.mix(FACETS.out.versions)

    //
    // Synchronization Point
    //
    ch_all_finished = ASCAT.out.segments
        .mix(FACETS.out.vcf, SEQUENZA_RUN.out.results)
        .collect()

    //
    // Check for ACESeq Output Directory
    //
    CHECK_ACESEQ_DIR(ch_all_finished)

    emit:
    versions
}