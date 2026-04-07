#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    VALIDATE INPUTS
*/

// Check input path parameters to see if they exist
def dbsnp_chr_prefix_rename_path = params.dbsnp_chr_prefix_rename ?: "${projectDir}/assets/dbsnp_chr_prefix_rename.tsv"
def checkPathParamList = [ 
    params.input, 
    params.aceseq_manifest,
    params.fasta, 
    params.fasta_fai, 
    params.ascat_alleles, 
    params.ascat_loci, 
    params.ascat_gc, 
    params.ascat_rt, 
    params.dbsnp, 
    params.dbsnp_tbi, 
    params.gcmapdir_frag,
    params.hapmap_sites,
    params.pon_dryclean,
    params.bwa_index,
    params.gridss_pon,
    params.brass_cache_dir,
    params.brass_depth,
    params.brass_viral,
    params.brass_gcbins,
    params.brass_cytoband,
    params.brass_centtel,
    params.brass_np,
    params.indel_mask,
    params.germ_sv_db,
    params.simple_seq_db,
    dbsnp_chr_prefix_rename_path
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.aceseq_manifest) { manifest_input = file(params.aceseq_manifest, checkIfExists: true) } else { exit 1, 'ACEseq manifest not specified!' }

/*
    CHANNEL SETUP
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
fasta                       = params.fasta                      ? Channel.fromPath(params.fasta).first()                        : Channel.empty()
fasta_fai                   = params.fasta_fai                  ? Channel.fromPath(params.fasta_fai).first()                  : Channel.empty()
target_regs                 = params.target_regs                ? Channel.fromPath(params.target_regs).collect()                : Channel.value([])
dbsnp                       = params.dbsnp                      ? Channel.fromPath(params.dbsnp).first()                        : Channel.empty()
dbsnp_tbi                   = params.dbsnp_tbi                  ? Channel.fromPath(params.dbsnp_tbi).first()                    : Channel.empty()
dbsnp_chr_prefix_rename     = Channel.fromPath(dbsnp_chr_prefix_rename_path).first()

// Ascat
ascat_alleles               = params.ascat_alleles              ? Channel.fromPath(params.ascat_alleles).collect()              : Channel.empty()
ascat_loci                  = params.ascat_loci                 ? Channel.fromPath(params.ascat_loci).collect()                 : Channel.empty()
ascat_gc                    = params.ascat_gc                   ? Channel.fromPath(params.ascat_gc).collect()                   : Channel.empty()
ascat_rt                    = params.ascat_rt                   ? Channel.fromPath(params.ascat_rt).collect()                   : Channel.empty()

// Facets
facets_annotation_bed       = params.facets_annotation_bed      ? Channel.fromPath(params.facets_annotation_bed).collect()      : Channel.value([])

// FragCounter
gcmapdir_frag               = params.gcmapdir_frag              ? Channel.fromPath(params.gcmapdir_frag).collect()              : Channel.empty()

// Delly
delly_blacklist             = params.delly_blacklist            ? Channel.fromPath(params.delly_blacklist).collect()            : Channel.value([])

// Gridss
bwa_index                   = params.bwa_index                  ? Channel.fromPath(params.bwa_index).collect()                  : Channel.empty()
gridss_blacklist            = params.gridss_blacklist           ? Channel.fromPath(params.gridss_blacklist).collect()           : Channel.empty()
gridss_pon                  = params.gridss_pon                 ? Channel.fromPath(params.gridss_pon).collect()                 : Channel.empty()

// Svaba
indel_mask                  = params.indel_mask                 ? Channel.fromPath(params.indel_mask).collect()                 : Channel.empty()
germ_sv_db                  = params.germ_sv_db                 ? Channel.fromPath(params.germ_sv_db).collect()                 : Channel.empty()
simple_seq_db               = params.simple_seq_db              ? Channel.fromPath(params.simple_seq_db).collect()              : Channel.empty()

// Brass
brass_cache_dir             = params.brass_cache_dir            ? Channel.fromPath(params.brass_cache_dir).collect()            : Channel.empty()
brass_depth                 = params.brass_depth                ? Channel.fromPath(params.brass_depth).collect()                : Channel.empty()
brass_viral                 = params.brass_viral                ? Channel.fromPath(params.brass_viral).collect()                : Channel.empty()
brass_microbe               = params.brass_microbe              ? Channel.fromPath(params.brass_microbe).collect()              : Channel.empty()
brass_gcbins                = params.brass_gcbins               ? Channel.fromPath(params.brass_gcbins).collect()               : Channel.empty()
brass_cytoband              = params.brass_cytoband             ? Channel.fromPath(params.brass_cytoband).collect()             : Channel.empty()
brass_centtel               = params.brass_centtel              ? Channel.fromPath(params.brass_centtel).collect()              : Channel.empty()
brass_np                    = params.brass_np                   ? Channel.fromPath(params.brass_np).collect()                   : Channel.empty()
brass_np_tbi                = params.brass_np_tbi               ? Channel.fromPath(params.brass_np_tbi).collect()               : Channel.empty()

// HetPileups
hapmap_sites                = params.hapmap_sites               ? Channel.fromPath(params.hapmap_sites).collect()               : Channel.empty()

// Dryclean
pon_dryclean                = params.pon_dryclean               ? Channel.fromPath(params.pon_dryclean).collect()               : Channel.empty()
blacklist_path_dryclean     = params.blacklist_path_dryclean    ? Channel.fromPath(params.blacklist_path_dryclean).collect()   : Channel.empty()
germline_file_dryclean      = params.germline_file_dryclean     ? Channel.fromPath(params.germline_file_dryclean).collect()    : Channel.empty()

/*
    IMPORT SUBWORKFLOWS
*/
include { INPUT_PREP } from '../subworkflows/input_prep.nf'
include { ZERO_SHOT_CNV_CALL } from '../subworkflows/zeroshot_cnv_calling.nf'
include { BAM_HETPILEUPS } from '../subworkflows/bam_hetpileups/main.nf'
include { COVERAGE_JABBA } from '../subworkflows/coverage_jabba.nf'
include { BREAKPOINT_ESTIMATOR } from '../subworkflows/breakpoint_estimation.nf'
include { BCFTOOLS_ANNOTATE } from '../modules/nf-core/bcftools/annotate/main.nf'
include { VALIDATE_ACESEQ_MANIFEST } from '../modules/local/check_aceseq/main.nf'

/*
    MAIN WORKFLOW
*/
workflow SOMASV_CALLER {
    versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    prepared_samples = INPUT_PREP(ch_input)

    VALIDATE_ACESEQ_MANIFEST(ch_input, manifest_input)
    manifest_ready = VALIDATE_ACESEQ_MANIFEST.out.ready
    samples = prepared_samples
        .combine(manifest_ready)
        .map { meta, control, control_index, tumor, tumor_index, _ ->
            [meta, control, control_index, tumor, tumor_index]
        }

    //println "The samples: "
    // ch_samples.view()

    //
    // MODULE: Normalize dbSNP contig names by adding chr-prefix before the downstream use
    //
    ch_dbsnp_chr_prefix = dbsnp
        .combine(dbsnp_tbi)
        .combine(dbsnp_chr_prefix_rename)
        .map { dbsnp_vcf, dbsnp_index, rename_map ->
            [[id: 'dbsnp'], dbsnp_vcf, dbsnp_index, [], [], [], [], rename_map]
        }

    BCFTOOLS_ANNOTATE(ch_dbsnp_chr_prefix)
    dbsnp_chr      = BCFTOOLS_ANNOTATE.out.vcf.map { meta, vcf -> vcf }.first()
    dbsnp_chr_tbi  = BCFTOOLS_ANNOTATE.out.tbi.map { meta, tbi -> tbi }.first()
    versions       = versions.mix(BCFTOOLS_ANNOTATE.out.versions_bcftools)

    //
    // SUBWORKFLOW: Run zero-shot CNV calling with ASCAT, SEQUENZA, and FACETS
    //  
    ZERO_SHOT_CNV_CALL(
        samples, 
        fasta, 
        fasta_fai,
        ascat_alleles, 
        ascat_loci, 
        ascat_gc, 
        ascat_rt, 
        target_regs,
        dbsnp_chr,
        dbsnp_chr_tbi,
        facets_annotation_bed
    )
    ascat_purityploidy  = ZERO_SHOT_CNV_CALL.out.ascat_purityploidy
    versions            = versions.mix(ZERO_SHOT_CNV_CALL.out.versions)

    //
    // SUBWORKFLOW: Generate GC- and mappability-corrected denoised genomic coverage files for JAbBA
    //
    COVERAGE_JABBA(
        samples, 
        params.midpoint_frag, 
        params.windowsize_frag, 
        gcmapdir_frag, 
        params.minmapq_frag, 
        params.paired_frag, 
        params.exome_frag, 
        pon_dryclean, 
        params.center_dryclean, 
        params.cbs_dryclean, 
        params.cnsignif_dryclean, 
        params.wholeGenome_dryclean, 
        params.blacklist_dryclean, 
        blacklist_path_dryclean, 
        params.germline_filter_dryclean, 
        germline_file_dryclean, 
        params.field_dryclean, 
        params.build_dryclean
    )
    versions    = versions.mix(COVERAGE_JABBA.out.versions)

    //
    // SUBWORKFLOW: Generate heterozygous pileups file for JabBa
    //
    BAM_HETPILEUPS(
        samples, 
        params.filter_hets, 
        params.max_depth_hets,
        hapmap_sites
    )
    versions = versions.mix(BAM_HETPILEUPS.out.versions)
    sites_from_het_pileups_wgs = Channel.empty().mix(BAM_HETPILEUPS.out.het_pileups_wgs)

    //
    // SUBWORKFLOW: Estimate genomic breakpoints for JabBa
    //
    BREAKPOINT_ESTIMATOR (
        samples,
        ascat_purityploidy,
        fasta,
        fasta_fai,
        dbsnp_chr,
        dbsnp_chr_tbi,
        bwa_index,
        indel_mask,
        germ_sv_db,
        simple_seq_db,
        params.svaba_error_rate,
        gridss_blacklist,
        gridss_pon,
        delly_blacklist,
        params.delly_mode,
        params.delly_min_svqual,
        params.delly_altaf,
        params.delly_min_svsize,
        params.delly_max_svsize,
        params.delly_ratio_geno,
        params.delly_pass,
        params.delly_tags,
        params.delly_cov,
        params.delly_ctrl_contamination,
        params.delly_geno_qual,
        params.delly_rddel,
        params.delly_rddup,
        brass_cache_dir,
        brass_depth,
        brass_viral,
        brass_microbe,
        brass_gcbins,
        brass_cytoband,
        brass_centtel,
        brass_np,
        brass_np_tbi,
        params.brass_sp,
        params.brass_protocol,
        params.brass_genome,
        params.brass_minreads,
        params.brass_mincn
    )
    versions = versions.mix(BREAKPOINT_ESTIMATOR.out.versions)
}
