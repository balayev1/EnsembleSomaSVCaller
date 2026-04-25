#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Define genomic breakpoints for JAbBA
//

// Delly
include { BAM_DELLY } from '../subworkflows/bam_delly/main.nf'
include { BCFTOOLS_VIEW } from '../modules/nf-core/bcftools/view/main.nf'

// Gridss
include { GRIDSS_SV_CALLING } from '../subworkflows/bam_gridss/main.nf'
include { GRIDSS_SOMATIC_FILTER_STEP } from '../subworkflows/bam_gridss/main.nf'

// Consensus merging and filtering
include { VCF_TO_PLAIN as DELLY_VCF_TO_PLAIN } from '../modules/local/vcf_to_plain/main.nf'
include { VCF_TO_PLAIN as GRIDSS_VCF_TO_PLAIN } from '../modules/local/vcf_to_plain/main.nf'
include { VCF_TO_PLAIN as BRASS_VCF_TO_PLAIN } from '../modules/local/vcf_to_plain/main.nf'
include { VCF_TO_PLAIN as SVABA_VCF_TO_PLAIN } from '../modules/local/vcf_to_plain/main.nf'
include { POST_SURVIVOR_PROCESSING as SURVIVOR_POST_FILTER_PROCESSING } from '../modules/local/post_survivor_processing/main.nf'
include { BCFTOOLS_REHEADER as SURVIVOR_BCFTOOLS_REHEADER } from '../modules/nf-core/bcftools/reheader/main.nf'
include { SURVIVOR_MERGE as BREAKPOINTS_SURVIVOR_MERGE } from '../modules/nf-core/survivor/merge/main.nf'
include { BEDTOOLS_SLOP as INDEL_MASK_SLOP } from '../modules/nf-core/bedtools/slop/main.nf'
include { BEDTOOLS_INTERSECT as FILTER_BLACKLIST_REGIONS } from '../modules/nf-core/bedtools/intersect/main.nf'

// Svaba
include { SVABA } from '../modules/local/svaba/main.nf'

// Brass
include { BAM_STATS as NORMAL_BAM_STATS } from '../modules/local/bam_stats/main.nf'
include { BAM_STATS as TUMOR_BAM_STATS }  from '../modules/local/bam_stats/main.nf'
include { BRASS }         from '../modules/local/brass/main.nf'
include { TRANSFORM_ASCAT_STATS } from '../modules/local/transform_ascat_stats/main.nf'

// Genomic breakpoint estimation workflow for JAbBA
workflow BREAKPOINT_ESTIMATOR {
    take:
    samples       // channel: [ val(meta), control, control_index, tumor, tumor_index ]
    ascat_purityploidy
    fasta            
    fai
    dbsnp
    dbsnp_tbi
    bwa_index
    indel_mask
    germ_sv_db
    simple_seq_db
    chrom_sizes
    svaba_error_rate
    survivor_max_distance_breakpoints
    survivor_min_supporting_callers
    survivor_account_for_type
    survivor_account_for_sv_strands
    survivor_estimate_distanced_by_sv_size
    survivor_min_sv_size
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
    brass_cache_dir
    brass_depth
    brass_viral
    brass_microbes
    brass_gcbins
    brass_cytoband
    brass_centtel
    brass_np
    brass_np_tbi
    brass_species
    brass_protocol
    brass_assembly
    brass_min_reads
    brass_mincn

    main:
    versions = Channel.empty()
    empty_intersect_sizes = Channel.value([[id: 'no_sizes'], []])
    simple_seq_bed = simple_seq_db.map { it instanceof List ? it[0] : it }
    indel_mask_bed = indel_mask.map { it instanceof List ? it[0] : it }

    //
    // Prepare reference inputs for DELLY
    //
    ch_fasta = fasta.map { it -> [ [id:'fasta'], it ] }
    ch_fai   = fai.map   { it -> [ [id:'fasta_fai'], it ] }

    //
    // Transform samples into DELLY input format
    //
    delly_input = samples
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
        samples,
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
    gridss_vcf_all  = GRIDSS_SOMATIC_FILTER_STEP.out.somatic_all
    versions        = versions.mix(GRIDSS_SOMATIC_FILTER_STEP.out.versions)

    //
    // Run SVABA
    //
    SVABA(
        samples,
        fasta,
        fai,
        bwa_index,
        dbsnp,
        dbsnp_tbi,
        indel_mask,
        germ_sv_db,
        simple_seq_db,
        svaba_error_rate
    )
    svaba_vcf_som_sv            = SVABA.out.som_sv              // high-confidence somatic sv
    svaba_vcf_unfiltered_som_sv = SVABA.out.unfiltered_som_sv   // unfiltered somatic sv
    versions                    = versions.mix(SVABA.out.versions)

    //
    // Transform ASCAT stats for BRASS
    //
    TRANSFORM_ASCAT_STATS(ascat_purityploidy)
    versions = versions.mix(TRANSFORM_ASCAT_STATS.out.versions)

    //
    // Run bam stats to make .bas file with BAM summary statistics
    //
    normal_bams         = samples.map { meta, control, control_index, tumor, tumor_index -> 
        [meta, control, control_index] 
    }

    tumor_bams          = samples.map { meta, control, control_index, tumor, tumor_index -> 
        [meta, tumor, tumor_index] 
    }

    NORMAL_BAM_STATS(normal_bams, fai)
    versions = versions.mix(NORMAL_BAM_STATS.out.versions)

    TUMOR_BAM_STATS(tumor_bams, fai)
    versions = versions.mix(TUMOR_BAM_STATS.out.versions)

    brass_input = samples
        .join(TUMOR_BAM_STATS.out.bas)
        .join(NORMAL_BAM_STATS.out.bas)
        .join(TRANSFORM_ASCAT_STATS.out.brass_stats)
        .map { meta, control, control_index, tumor, tumor_index, 
               t_bas, n_bas, ascat_sum ->
               [ meta, control, control_index, tumor, tumor_index, t_bas, n_bas, ascat_sum ]
        }

    //
    // Run BRASS
    //
    BRASS(
        brass_input,
        fasta,
        fai,
        brass_cache_dir,
        brass_depth,
        brass_viral,
        brass_microbes,
        brass_gcbins,
        brass_cytoband,
        brass_centtel,
        brass_np,
        brass_np_tbi,
        brass_species,
        brass_protocol,
        brass_assembly,
        brass_min_reads,
        brass_mincn
    )
    brass_bedpe = BRASS.out.bedpe
    brass_vcf   = BRASS.out.vcf
    versions    = versions.mix(BRASS.out.versions)

    //
    // Convert compressed caller VCFs to plain-text for SURVIVOR
    //
    DELLY_VCF_TO_PLAIN(delly_vcf)
    GRIDSS_VCF_TO_PLAIN(gridss_vcf_hc)
    BRASS_VCF_TO_PLAIN(brass_vcf)
    SVABA_VCF_TO_PLAIN(svaba_vcf_som_sv)
    versions = versions
        .mix(DELLY_VCF_TO_PLAIN.out.versions_gzip)
        .mix(GRIDSS_VCF_TO_PLAIN.out.versions_gzip)
        .mix(BRASS_VCF_TO_PLAIN.out.versions_gzip)
        .mix(SVABA_VCF_TO_PLAIN.out.versions_gzip)

    survivor_merge_input = DELLY_VCF_TO_PLAIN.out.vcf
        .join(GRIDSS_VCF_TO_PLAIN.out.vcf)
        .join(BRASS_VCF_TO_PLAIN.out.vcf)
        .join(SVABA_VCF_TO_PLAIN.out.vcf)
        .map { meta, delly_plain_vcf, gridss_plain_vcf, brass_plain_vcf, svaba_plain_vcf ->
            [meta, [delly_plain_vcf, gridss_plain_vcf, brass_plain_vcf, svaba_plain_vcf]]
        }

    //
    // Merge consensus breakpoints across callers
    //
    BREAKPOINTS_SURVIVOR_MERGE(
        survivor_merge_input,
        survivor_max_distance_breakpoints,
        survivor_min_supporting_callers,
        survivor_account_for_type,
        survivor_account_for_sv_strands,
        survivor_estimate_distanced_by_sv_size,
        survivor_min_sv_size
    )
    survivor_consensus_vcf = BREAKPOINTS_SURVIVOR_MERGE.out.vcf
    versions               = versions.mix(BREAKPOINTS_SURVIVOR_MERGE.out.versions_survivor)

    //
    // Build padded blacklist regions per sample
    //
    INDEL_MASK_SLOP(
        survivor_consensus_vcf
            .combine(indel_mask_bed)
            .map { meta, consensus_vcf, indel_mask_bed_file -> [meta, indel_mask_bed_file] },
        chrom_sizes
    )
    versions = versions.mix(INDEL_MASK_SLOP.out.versions_bedtools)

    //
    // Filter consensus breakpoints against blacklist regions
    //
    FILTER_BLACKLIST_REGIONS(
        survivor_consensus_vcf
            .join(INDEL_MASK_SLOP.out.bed)
            .map { meta, consensus_vcf_file, blacklist_bed -> [meta, consensus_vcf_file, blacklist_bed] },
        empty_intersect_sizes
    )
    survivor_blacklist_filtered_vcf = FILTER_BLACKLIST_REGIONS.out.intersect
    versions                        = versions.mix(FILTER_BLACKLIST_REGIONS.out.versions_bedtools)

    //
    // Remove chrY breakpoints from female samples after blacklist filtering
    //
    SURVIVOR_POST_FILTER_PROCESSING(survivor_blacklist_filtered_vcf)
    survivor_postprocessed_vcf = SURVIVOR_POST_FILTER_PROCESSING.out.vcf
    versions                   = versions.mix(SURVIVOR_POST_FILTER_PROCESSING.out.versions)

    //
    // Restore contig lengths from the reference .fai before JaBbA seqname coercion
    //
    SURVIVOR_BCFTOOLS_REHEADER(
        survivor_postprocessed_vcf.map { meta, vcf -> [meta, vcf, [], []] },
        ch_fai
    )
    survivor_filtered_vcf = SURVIVOR_BCFTOOLS_REHEADER.out.vcf
    versions              = versions.mix(SURVIVOR_BCFTOOLS_REHEADER.out.versions_bcftools)

    emit:
    delly_filter_bcf
    delly_vcf
    gridss_vcf_hc
    gridss_vcf_all
    svaba_vcf_som_sv
    svaba_vcf_unfiltered_som_sv
    brass_bedpe
    brass_vcf
    survivor_consensus_vcf
    survivor_filtered_vcf
    versions
}
