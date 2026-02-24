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
    svaba_error_rate
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
    brass_genome_cache
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
        brass_genome_cache,
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

    emit:
    delly_filter_bcf
    delly_vcf
    gridss_vcf_hc
    gridss_vcf_all
    svaba_vcf_som_sv
    svaba_vcf_unfiltered_som_sv
    brass_bedpe
    brass_vcf
    versions
}