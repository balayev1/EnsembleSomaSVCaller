#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Calling CNVs using a zero-shot approach with FOUR tools ASCAT, SEQUENZA, FACETS, ACESEQ
//

include { ASCAT }                   from '../modules/nf-core/ascat/main.nf'
include { FACETS }                  from '../modules/local/facets/main.nf'
include { SEQUENZAUTILS_GCWIGGLE }  from '../modules/local/sequenza/main.nf'
include { SEQUENZA_PREP }           from '../modules/local/sequenza/main.nf'
include { SEQUENZA_MERGE }          from '../modules/local/sequenza/main.nf'
include { SEQUENZA_RUN }            from '../modules/local/sequenza/main.nf'
include { CHECK_ACESEQ_DIR }        from '../modules/local/check_aceseq/main.nf'

// Zero-shot CNV Calling Workflow
workflow ZERO_SHOT_CNV_CALL {
    take:
        ch_samples    // channel: [val(meta), [control], [control_index],[ tumor], [tumor_index]]
        fasta         // channel: [path(fasta)]
        fai             // channel: [path(fai)]
        ascat_alleles    // channel: [path(allele_res)]
        ascat_loci      // channel: [path(ascat_loci)]
        ascat_gc       // channel: [path(ascat_gc)]
        ascat_rt       // channel: [path(ascat_rt)]
        target_regs    // channel: [path(target_regs)]
        dbsnp          // channel: [path(dbsnp)]
        dbsnp_tbi      // channel: [path(dbsnp_tbi)]
        facets_annotation // channel: [path(facets_annotation_bed)]

    main:
        versions = Channel.empty()

        //
        // Run ASCAT
        //
        ASCAT(
            ch_samples,
            ascat_alleles,
            ascat_loci,
            target_regs,
            fasta,
            ascat_gc,
            ascat_rt
        )
        ascat_purityploidy  = ASCAT.out.purityploidy
        versions            = versions.mix(ASCAT.out.versions)

        //
        // Run Sequenzautils to generate GC Wiggle Reference
        //
        ch_gc_wiggle_input = Channel.value([id:'reference_genome']).combine(fasta)
        SEQUENZAUTILS_GCWIGGLE(
            ch_gc_wiggle_input,
            params.windowsize_sequtils
        )
        ch_wiggle_file = SEQUENZAUTILS_GCWIGGLE.out.wig.map { meta, wig -> wig }.first()
        versions = versions.mix(SEQUENZAUTILS_GCWIGGLE.out.versions)

        //
        // Run Sequenzautils to generate binned seqz files
        //
        ch_chromosomes = fai
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .filter { it =~ /^(chr)?([0-9]+|[XY])$/ } // Filters for main contigs only
        ch_sequenza_prep_input = ch_samples.combine(ch_chromosomes)
        SEQUENZA_PREP(
            ch_sequenza_prep_input,
            fasta,
            ch_wiggle_file,
            params.hom_sequtils,
            params.het_sequtils,
            params.het_f_sequtils,
            params.qlimit_sequtils,
            params.qformat_sequtils,
            params.rd_thr_sequtils,
            params.seqz_bin_size_sequtils
        )

        ch_to_merge = SEQUENZA_PREP.out.seqz_file.groupTuple(by: 0)
        SEQUENZA_MERGE(ch_to_merge)
        versions = versions.mix(SEQUENZA_MERGE.out.versions)
        
        // 
        // Run SEQUENZA
        //
        SEQUENZA_RUN( 
            SEQUENZA_MERGE.out.merged_seqz,
            params.sequenza_genome,
            params.sequenza_rd_thr_tumor,
            params.sequenza_rd_thr_normal,
            params.sequenza_purity_range,
            params.sequenza_ploidy_range
        )
        versions = versions.mix(SEQUENZA_RUN.out.versions)

        // 
        // Run FACETS
        //
        FACETS(
            ch_samples,
            dbsnp,
            dbsnp_tbi,
            target_regs,
            facets_annotation,
            params.facets_mapq,
            params.facets_baq,
            params.facets_mindepth,
            params.facets_maxdepth,
            params.facets_cval_pre,
            params.facets_cval_proc,
            params.facets_nbhd_snp,
            params.facets_genome,
            params.facets_count_orphans,
            params.facets_unmatched,
            params.facets_no_cov_plot
        )
        versions = versions.mix(FACETS.out.versions)

        //
        // Synchronization Point
        //
        ch_all_finished = ASCAT.out.segments
            .mix(FACETS.out.vcf, SEQUENZA_RUN.out.segments)
            .collect()

        //
        // Check for ACESeq Output Directory
        //
        CHECK_ACESEQ_DIR(ch_all_finished)

    emit:
        ascat_purityploidy
        versions
}