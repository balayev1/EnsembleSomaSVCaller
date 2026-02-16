#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Calling CNVs using a zero-shot approach with FOUR tools ASCAT, SEQUENZA, FACETS, ACESEQ
//

params.options = [:]

include { ASCAT }                 from './modules/nf-core/ascat/main.nf' addParams( options: params.options )
include { FACETS }                from './modules/local/facets.nf' addParams( options: params.options )
include { SEQUENZAUTILS_GCWIGGLE } from './modules/local/sequenza.nf' addParams( options: params.options )
include { SEQUENZA_RUN }          from './modules/local/sequenza.nf' addParams( options: params.options )
include { CHECK_ACESEQ_DIR }      from '../modules/local/check_aceseq.nf' addParams( options: params.options )

// Zero-shot CNV Calling Workflow
workflow ZERO_SHOT_CNV_CALL {
    take
        ch_samples,    // channel: [val(meta), tumor,tumor_bai, control, control_bai]
        fasta,         // channel: [path(fasta)]
        ascat_alleles,    // channel: [path(allele_res)]
        ascat_loci,      // channel: [path(ascat_loci)]
        ascat_gc,       // channel: [path(ascat_gc)]
        ascat_rt,       // channel: [path(ascat_rt)]
        target_regs,       // channel: [path(target_regs)]
        dbsnp, // channel: [path(dbsnp)]
        dbsnp_tbi, // channel: [path(dbsnp_tbi)]
        facets_annotation // channel: [path(facets_annotation_bed)]

    main:
    versions = Channel.empty()

    //
    // Run ASCAT
    //
    ASCAT (
        ch_samples,
        ascat_alleles,
        ascat_loci,
        target_regs,
        fasta,
        ascat_gc,
        ascat_rt
    )
    versions = versions.mix(ASCAT.out.versions)

    //
    // Run Sequenzautils to generate GC Wiggle Reference
    //
    ch_gc_wiggle_input = Channel.of( [ [id:'reference_genome'], fasta ] )
    SEQUENZAUTILS_GCWIGGLE(ch_gc_wiggle_input)
    ch_wiggle_file = SEQUENZAUTILS_GCWIGGLE.out.wig.map { meta, wig -> wig }
    versions = versions.mix(SEQUENZAUTILS_GCWIGGLE.out.versions)

    // 
    // Run SEQUENZA
    //
    SEQUENZA_RUN ( 
        ch_samples, 
        fasta, 
        ch_wiggle_file 
    )
    versions = versions.mix(SEQUENZA_RUN.out.versions)

    // 
    // Run FACETS
    //
    FACETS (
        ch_samples,
        dbsnp,
        dbsnp_tbi,
        target_regs,
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